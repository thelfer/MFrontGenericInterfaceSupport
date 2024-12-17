/*!
 * \file   Integrate.cxx
 * \brief
 * \author Thomas Helfer
 * \date   21/08/2018
 * \copyright Copyright (C) 2006-2018 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence
 * or the CECILL-A licence. A copy of thoses licences are delivered
 * with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include <map>
#include <tuple>
#include <thread>
#include <memory>
#include <cstdlib>
#include <cinttypes>
#include "MGIS/Raise.hxx"
#include "MGIS/ThreadPool.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

namespace mgis::behaviour::internals {

  static void allocate(MaterialDataManager& m,
                       const BehaviourIntegrationOptions& opts) {
    if (opts.integration_type !=
        IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR) {
      m.allocateArrayOfTangentOperatorBlocks();
    }
    if (opts.compute_speed_of_sound) {
      m.allocateArrayOfSpeedOfSounds();
    }
  }  // end of allocate

  static mgis::real encodeBehaviourIntegrationOptions(
      const BehaviourIntegrationOptions& opts) {
    if (opts.compute_speed_of_sound) {
      return 100 + static_cast<int>(opts.integration_type);
    }
    return static_cast<int>(opts.integration_type);
  }  // end of encodeBehaviourIntegrationOptions

  //! \brief a simple alias
  using Evaluator = std::tuple<size_type,  // offset
                               size_type,  // variable size
                               const real*>;

  /*
   * \brief uniform values are treated immediatly. For spatially variable
   * fields, we return the information needed to evaluate them
   */
  static std::vector<Evaluator> buildEvaluators(
      std::vector<real>& v,
      std::map<std::string, MaterialStateManager::FieldHolder>& values,
      const MaterialDataManager& m,
      const std::vector<Variable>& ds) {
    mgis::raise_if(v.size() != getArraySize(ds, m.b.hypothesis),
                   "buildEvaluators: ill allocated memory");
    // evaluators
    std::vector<Evaluator> evaluators;
    auto offset = mgis::size_type{};
    for (const auto& d : ds) {
      const auto s = getVariableSize(d, m.b.hypothesis);
      auto p = values.find(d.name);
      if (p == values.end()) {
        auto msg = std::string{"buildEvaluators: no variable named '" + d.name +
                               "' declared"};
        if (!values.empty()) {
          msg += "\nThe following variables were declared: ";
          for (const auto& variable : values) {
            msg += "\n- " + variable.first;
          }
        } else {
          msg += "\nNo variable declared.";
        }
        mgis::raise(msg);
      }
      auto set_values = [&](const auto& variable_values) {
        if (variable_values.size() == s) {
          // uniform value
          std::copy(variable_values.begin(), variable_values.end(),
                    v.begin() + offset);
        } else {
          evaluators.push_back(
              std::make_tuple(offset, s, variable_values.data()));
        }
      };
      if (std::holds_alternative<real>(p->second)) {
        if (d.type != Variable::SCALAR) {
          mgis::raise(
              "buildEvaluators: invalid type for "
              "variable '" +
              d.name + "'");
        }
        v[offset] = std::get<real>(p->second);
      } else if (std::holds_alternative<std::span<real>>(p->second)) {
        set_values(std::get<std::span<real>>(p->second));
      } else {
        set_values(std::get<std::vector<real>>(p->second));
      }
      offset += s;
    }
    return evaluators;
  }  // end of buildEvaluators

  static std::optional<Evaluator> buildMassDensityEvaluator(
      mgis::real& rho, MaterialStateManager& s) {
    if (!isMassDensityDefined(s)) {
      rho = 0;
      return {};
    }
    if (isMassDensityUniform(s)) {
      rho = std::get<mgis::real>(*(s.mass_density));
      return {};
    }
    auto check = [&s](const auto& values) {
      if (values.size() != s.n) {
        mgis::raise(
            "buildMassDensityEvaluator: invalid size for the arrray of the "
            "mass density (" +
            std::to_string(values.size()) + " values given for '" +
            std::to_string(s.n) + "'integration points)");
      }
    };
    if (std::holds_alternative<std::vector<mgis::real>>(*(s.mass_density))) {
      const auto& values = std::get<std::vector<mgis::real>>(*(s.mass_density));
      check(values);
      return std::make_tuple(0, 1, values.data());
    }
    const auto& values = std::get<std::span<mgis::real>>(*(s.mass_density));
    check(values);
    return std::make_tuple(0, 1, values.data());
  }  // end of buildMassDensityEvaluator

  static inline void applyEvaluators(std::vector<real>& values,
                                     const std::vector<Evaluator>& evs,
                                     const size_type i) {
    for (const auto& ev : evs) {
      const auto o = std::get<0>(ev);
      const auto s = std::get<1>(ev);
      if (s == 1) {
        const auto* const v = std::get<2>(ev) + i;
        values[o] = *v;
      } else {
        const auto* const v = std::get<2>(ev) + i * s;
        std::copy(v, v + s, values.begin() + o);
      }
    }
  }  // end of applyEvaluators

  /*!
   * \brief small structure in charge of storing evaluators associated with a
   * behaviour.
   */
  struct BehaviourEvaluators {
    /*!
     * \brief evaluators associated with the material properties at the
     * beginning of the time step
     */
    std::vector<Evaluator> mps0;
    /*!
     * \brief evaluators associated with the material properties at the
     * end of the time step
     */
    std::vector<Evaluator> mps1;
    /*!
     * \brief evaluators associated with the external state variables at the
     * beginning of the time step
     */
    std::vector<Evaluator> esvs0;
    /*!
     * \brief evaluators associated with the external state variables at the
     * end of the time step
     */
    std::vector<Evaluator> esvs1;
    /*!
     * \brief evaluator associated with the mass density at the beginning of
     * the time step
     */
    std::optional<Evaluator> rho0;
    /*!
     * \brief evaluator associated with the mass density at the end of
     * the time step
     */
    std::optional<Evaluator> rho1;
  };  // end of struct BehaviourEvaluators

  static inline BehaviourEvaluators buildBehaviourEvaluators(
      BehaviourIntegrationWorkSpace& ws, MaterialDataManager& m) {
    //
    auto evaluators = BehaviourEvaluators{};
    // treating uniform values
    evaluators.mps0 = internals::buildEvaluators(
        ws.mps0, m.s0.material_properties, m, m.b.mps);
    evaluators.mps1 = internals::buildEvaluators(
        ws.mps1, m.s1.material_properties, m, m.b.mps);
    evaluators.esvs0 = internals::buildEvaluators(
        ws.esvs0, m.s0.external_state_variables, m, m.b.esvs);
    evaluators.esvs1 = internals::buildEvaluators(
        ws.esvs1, m.s1.external_state_variables, m, m.b.esvs);
    evaluators.rho0 = internals::buildMassDensityEvaluator(ws.rho0, m.s0);
    evaluators.rho1 = internals::buildMassDensityEvaluator(ws.rho1, m.s1);
    return evaluators;
  }  // end of buildBehaviourEvaluator

  static inline void evaluate(
      mgis::behaviour::BehaviourIntegrationWorkSpace& ws,
      const BehaviourEvaluators& evaluators,
      const size_type i) {
    applyEvaluators(ws.mps0, evaluators.mps0, i);
    applyEvaluators(ws.mps1, evaluators.mps1, i);
    applyEvaluators(ws.esvs0, evaluators.esvs0, i);
    applyEvaluators(ws.esvs1, evaluators.esvs1, i);
    if (evaluators.rho0.has_value()) {
      ws.rho0 = *(std::get<2>(*(evaluators.rho0)) + i);
    }
    if (evaluators.rho1.has_value()) {
      ws.rho1 = *(std::get<2>(*(evaluators.rho1)) + i);
    }
  }

  static inline mgis::behaviour::BehaviourDataView initializeBehaviourDataView(
      mgis::behaviour::BehaviourIntegrationWorkSpace& ws) {
    auto v = mgis::behaviour::BehaviourDataView{};
    v.error_message = ws.error_message.data();
    v.s0.material_properties = ws.mps0.data();
    v.s1.material_properties = ws.mps1.data();
    v.s0.external_state_variables = ws.esvs0.data();
    v.s1.external_state_variables = ws.esvs1.data();
    v.s0.mass_density = &(ws.rho0);
    v.s1.mass_density = &(ws.rho1);
    v.s0.stored_energy = nullptr;
    v.s1.stored_energy = nullptr;
    v.s0.dissipated_energy = nullptr;
    v.s1.dissipated_energy = nullptr;
    return v;
  }  // end of initializeBehaviourDataView

  static inline void updateView(mgis::behaviour::BehaviourDataView& v,
                                const mgis::behaviour::MaterialDataManager& m,
                                const size_type i) {
    // strides
    const auto g_stride = m.s0.gradients_stride;
    const auto t_stride = m.s0.thermodynamic_forces_stride;
    const auto isvs_stride = m.s0.internal_state_variables_stride;
    const auto computes_stored_energy = m.b.computesStoredEnergy;
    const auto computes_dissipated_energy = m.b.computesDissipatedEnergy;
    v.speed_of_sound = m.speed_of_sound.data() + i;
    v.s0.gradients = m.s0.gradients.data() + g_stride * i;
    v.s1.gradients = m.s1.gradients.data() + g_stride * i;
    v.s0.thermodynamic_forces = m.s0.thermodynamic_forces.data() + t_stride * i;
    v.s1.thermodynamic_forces = m.s1.thermodynamic_forces.data() + t_stride * i;
    v.s0.internal_state_variables =
        m.s0.internal_state_variables.data() + isvs_stride * i;
    v.s1.internal_state_variables =
        m.s1.internal_state_variables.data() + isvs_stride * i;
    if (computes_stored_energy) {
      v.s0.stored_energy = m.s0.stored_energies.data() + i;
      v.s1.stored_energy = m.s1.stored_energies.data() + i;
    }
    if (computes_dissipated_energy) {
      v.s0.dissipated_energy = m.s0.dissipated_energies.data() + i;
      v.s1.dissipated_energy = m.s1.dissipated_energies.data() + i;
    }
  }  // end of updateView

  static inline void checkIntegrationPointsRange(
      const mgis::behaviour::MaterialDataManager& m,
      const size_type b,
      const size_type e) {
    if (b > e) {
      mgis::raise(
          "checkIntegrationPointsRange: "
          "invalid range (lower bound is greated than the upper bound)");
    }
    if (e > m.n) {
      mgis::raise(
          "checkIntegrationPointsRange: "
          "invalid upper bound ('" +
          std::to_string(e) + "')");
    }
  }  // end of checkIntegrationPointsRange

  /*!
   * \brief execute the given initialize function over a range of integration
   * points.
   */
  static BehaviourIntegrationResult executeInitializeFunction(
      MaterialDataManager& m,
      const BehaviourInitializeFunction p,
      const mgis::size_type b,
      const mgis::size_type e) {
    // workspace
    auto& ws = m.getBehaviourIntegrationWorkSpace();
    auto v = internals::initializeBehaviourDataView(ws);
    auto behaviour_evaluators = internals::buildBehaviourEvaluators(ws, m);
    v.rdt = nullptr;
    // loop over integration points
    auto r = BehaviourIntegrationResult{};
    for (auto i = b; i != e; ++i) {
      internals::evaluate(ws, behaviour_evaluators, i);
      internals::updateView(v, m, i);
      v.dt = mgis::real{};
      const auto ri = (p.f)(&v, nullptr);
      if (ri != 0) {
        r.n = i;
        v.error_message[511] = '\0';
        r.error_message = std::string(v.error_message);
        return r;
      }
    }
    return r;
  }  // end of executeInitializeFunction

  /*!
   * \brief execute the given initialize function over a range of integration
   * points.
   */
  static BehaviourIntegrationResult executeInitializeFunction(
      MaterialDataManager& m,
      const BehaviourInitializeFunction p,
      std::span<const real> inputs,
      const mgis::size_type inputs_stride,
      const mgis::size_type b,
      const mgis::size_type e) {
    // workspace
    auto& ws = m.getBehaviourIntegrationWorkSpace();
    auto v = internals::initializeBehaviourDataView(ws);
    auto behaviour_evaluators = internals::buildBehaviourEvaluators(ws, m);
    v.rdt = nullptr;
    // loop over integration points
    auto r = BehaviourIntegrationResult{};
    const auto* const inputs_values = inputs.data();
    for (auto i = b; i != e; ++i) {
      internals::evaluate(ws, behaviour_evaluators, i);
      internals::updateView(v, m, i);
      v.dt = mgis::real{};
      const auto ri = (p.f)(&v, inputs_values + inputs_stride * i);
      if (ri != 0) {
        r.n = i;
        v.error_message[511] = '\0';
        r.error_message = std::string(v.error_message);
        return r;
      }
    }
    return r;
  }  // end of executeInitializeFunction

  /*!
   * \brief perform the integration of the behaviour over a range of integration
   * points.
   */
  static BehaviourIntegrationResult integrate(
      MaterialDataManager& m,
      const BehaviourIntegrationOptions& opts,
      const real dt,
      const size_type b,
      const size_type e) {
    // workspace
    auto& ws = m.getBehaviourIntegrationWorkSpace();
    auto v = internals::initializeBehaviourDataView(ws);
    auto behaviour_evaluators = internals::buildBehaviourEvaluators(ws, m);
    // loop over integration points
    auto r = BehaviourIntegrationResult{};
    auto rdt0 = r.time_step_increase_factor;
    const real Ke = encodeBehaviourIntegrationOptions(opts);
    real bopts[Behaviour::nopts + 1];  // option passed to the behaviour
    for (auto i = b; i != e; ++i) {
      internals::evaluate(ws, behaviour_evaluators, i);
      internals::updateView(v, m, i);
      auto rdt = rdt0;
      v.error_message[0] = '\0';
      v.rdt = &rdt;
      v.dt = dt;
      if ((opts.integration_type !=
           IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR) &&
          (m.K_stride != 0)) {
        v.K = m.K.data() + m.K_stride * i;
      } else {
        v.K = &bopts[0];
      }
      v.K[0] = Ke;
      const auto ri = integrate(v, m.b);
      r.exit_status = std::min(ri, r.exit_status);
      r.time_step_increase_factor = std::min(rdt, r.time_step_increase_factor);
      if (ri == 0) {
        r.n = i;
      } else if (ri == -1) {
        r.n = i;
        v.error_message[511] = '\0';
        r.error_message = std::string(v.error_message);
        return r;
      }
    }
    return r;
  }  // end of integrate

  /*!
   * \brief execute the given post-processing over a range of integration
   * points.
   */
  static BehaviourIntegrationResult executePostProcessing(
      std::span<real> outputs,
      MaterialDataManager& m,
      const BehaviourPostProcessing p,
      const mgis::size_type outputs_stride,
      const mgis::size_type b,
      const mgis::size_type e) {
    // workspace
    auto& ws = m.getBehaviourIntegrationWorkSpace();
    auto v = internals::initializeBehaviourDataView(ws);
    auto behaviour_evaluators = internals::buildBehaviourEvaluators(ws, m);
    v.rdt = nullptr;
    // loop over integration points
    auto r = BehaviourIntegrationResult{};
    auto* const outputs_values = outputs.data();
    for (auto i = b; i != e; ++i) {
      internals::evaluate(ws, behaviour_evaluators, i);
      internals::updateView(v, m, i);
      v.dt = mgis::real{};
      const auto ri = (p.f)(outputs_values + outputs_stride * i, &v);
      if (ri != 0) {
        r.n = i;
        v.error_message[511] = '\0';
        r.error_message = std::string(v.error_message);
        return r;
      }
    }
    return r;
  }  // end of executePostProcessing

}  // namespace mgis::behaviour::internals

namespace mgis::behaviour {

  BehaviourIntegrationResult::BehaviourIntegrationResult() = default;

  BehaviourIntegrationResult::BehaviourIntegrationResult(
      BehaviourIntegrationResult&&) = default;

  BehaviourIntegrationResult::BehaviourIntegrationResult(
      const BehaviourIntegrationResult&) = default;

  BehaviourIntegrationResult& BehaviourIntegrationResult::operator=(
      BehaviourIntegrationResult&&) = default;

  BehaviourIntegrationResult& BehaviourIntegrationResult::operator=(
      const BehaviourIntegrationResult&) = default;

  BehaviourIntegrationResult::~BehaviourIntegrationResult() = default;

  MultiThreadedBehaviourIntegrationResult::
      MultiThreadedBehaviourIntegrationResult() = default;

  MultiThreadedBehaviourIntegrationResult::
      MultiThreadedBehaviourIntegrationResult(
          MultiThreadedBehaviourIntegrationResult&&) = default;

  MultiThreadedBehaviourIntegrationResult::
      MultiThreadedBehaviourIntegrationResult(
          const MultiThreadedBehaviourIntegrationResult&) = default;

  MultiThreadedBehaviourIntegrationResult&
  MultiThreadedBehaviourIntegrationResult::operator=(
      MultiThreadedBehaviourIntegrationResult&&) = default;

  MultiThreadedBehaviourIntegrationResult&
  MultiThreadedBehaviourIntegrationResult::operator=(
      const MultiThreadedBehaviourIntegrationResult&) = default;

  MultiThreadedBehaviourIntegrationResult::
      ~MultiThreadedBehaviourIntegrationResult() = default;

  static const BehaviourInitializeFunction& getBehaviourInitializeFunction(
      const Behaviour& b, const std::string_view n) {
    const auto p = b.initialize_functions.find(n);
    if (p == b.initialize_functions.end()) {
      mgis::raise(
          "executeInitializeFunction: "
          "no initialize function named '" +
          std::string{n} + "'");
    }
    return p->second;
  }  // end of getBehaviourInitializeFunction

  int executeInitializeFunction(BehaviourDataView& d,
                                const Behaviour& b,
                                const std::string_view n) {
    const auto ifct = getBehaviourInitializeFunction(b, n);
    if (!ifct.inputs.empty()) {
      mgis::raise(
          "executeInitializeFunction: "
          "invalid size of the inputs '" +
          std::string{n} + "'");
    }
    const auto& f = ifct.f;
    return f(&d, nullptr);
  }  // end of executeInitializeFunction

  int executeInitializeFunction(BehaviourDataView& d,
                                const Behaviour& b,
                                const std::string_view n,
                                std::span<const real> inputs) {
    const auto ifct = getBehaviourInitializeFunction(b, n);
    if (inputs.size() != getArraySize(ifct.inputs, b.hypothesis)) {
      mgis::raise(
          "executeInitializeFunction: "
          "invalid size of the inputs '" +
          std::string{n} + "'");
    }
    const auto& f = ifct.f;
    return f(&d, inputs.data());
  }  // end of executeInitializeFunction

  BehaviourIntegrationResult executeInitializeFunction(MaterialDataManager& m,
                                                       const std::string_view n,
                                                       const size_type b,
                                                       const size_type e) {
    internals::checkIntegrationPointsRange(m, b, e);
    const auto& ifct = getBehaviourInitializeFunction(m.b, n);
    if (!ifct.inputs.empty()) {
      mgis::raise(
          "executeInitializeFunction: "
          "invalid size of the inputs '" +
          std::string{n} + "'");
    }
    return internals::executeInitializeFunction(m, ifct, b, e);
  }  // end of executeInitializeFunction

  BehaviourIntegrationResult executeInitializeFunction(
      MaterialDataManager& m,
      const std::string_view n,
      std::span<const real> inputs,
      const size_type b,
      const size_type e) {
    const auto& ifct = getBehaviourInitializeFunction(m.b, n);
    const auto istride = getArraySize(ifct.inputs, m.b.hypothesis);
    internals::checkIntegrationPointsRange(m, b, e);
    if ((inputs.size() != m.n * istride) && (inputs.size() != istride)) {
      mgis::raise(
          "executeInitializeFunction: "
          "invalid size of the inputs '" +
          std::string{n} + "'");
    }
    if (inputs.size() == istride) {
      return internals::executeInitializeFunction(m, ifct, inputs, 0, b, e);
    }
    return internals::executeInitializeFunction(m, ifct, inputs, istride, b, e);
  }  // end of executeInitializeFunction

  BehaviourIntegrationResult executeInitializeFunction(
      MaterialDataManager& m, const std::string_view n) {
    return executeInitializeFunction(m, n, 0, m.n);
  }  // end of executeInitializeFunction

  BehaviourIntegrationResult executeInitializeFunction(
      MaterialDataManager& m,
      const std::string_view n,
      std::span<const real> inputs) {
    return executeInitializeFunction(m, n, inputs, 0, m.n);
  }  // end of executeInitializeFunction

  MultiThreadedBehaviourIntegrationResult executeInitializeFunction(
      mgis::ThreadPool& p, MaterialDataManager& m, const std::string_view n) {
    const auto& ifct = getBehaviourInitializeFunction(m.b, n);
    if (!ifct.inputs.empty()) {
      mgis::raise(
          "executeInitializeFunction: "
          "invalid size of the inputs '" +
          std::string{n} + "'");
    }
    m.setThreadSafe(true);
    // get number of threads
    const auto nth = p.getNumberOfThreads();
    const auto d = m.n / nth;
    const auto r = m.n % nth;
    size_type b = 0;
    std::vector<std::future<ThreadedTaskResult<BehaviourIntegrationResult>>>
        tasks;
    tasks.reserve(nth);
    for (size_type i = 0; i != r; ++i) {
      tasks.push_back(p.addTask([&m, &ifct, b, d] {
        return internals::executeInitializeFunction(m, ifct, b, b + d + 1);
      }));
      b += d + 1;
    }
    for (size_type i = r; i != nth; ++i) {
      tasks.push_back(p.addTask([&m, &ifct, b, d] {
        return internals::executeInitializeFunction(m, ifct, b, b + d);
      }));
      b += d;
    }
    auto res = MultiThreadedBehaviourIntegrationResult{};
    for (auto& t : tasks) {
      const auto& ri = *(t.get());
      res.exit_status = std::min(res.exit_status, ri.exit_status);
      res.results.push_back(ri);
    }
    return res;
  }  // end of executeInitializeFunction

  MultiThreadedBehaviourIntegrationResult executeInitializeFunction(
      mgis::ThreadPool& p,
      MaterialDataManager& m,
      const std::string_view n,
      std::span<const real> inputs) {
    const auto& ifct = getBehaviourInitializeFunction(m.b, n);
    const auto istride = getArraySize(ifct.inputs, m.b.hypothesis);
    if ((inputs.size() != m.n * istride) && (inputs.size() != istride)) {
      mgis::raise(
          "executeInitializeFunction: "
          "invalid size of the inputs '" +
          std::string{n} + "'");
    }
    // effective stride
    const auto estride = (inputs.size() == istride) ? 0 : istride;
    m.setThreadSafe(true);
    // get number of threads
    const auto nth = p.getNumberOfThreads();
    const auto d = m.n / nth;
    const auto r = m.n % nth;
    size_type b = 0;
    std::vector<std::future<ThreadedTaskResult<BehaviourIntegrationResult>>>
        tasks;
    tasks.reserve(nth);
    for (size_type i = 0; i != r; ++i) {
      tasks.push_back(p.addTask([&inputs, &m, &ifct, estride, b, d] {
        return internals::executeInitializeFunction(m, ifct, inputs, estride, b,
                                                    b + d + 1);
      }));
      b += d + 1;
    }
    for (size_type i = r; i != nth; ++i) {
      tasks.push_back(p.addTask([&inputs, &m, &ifct, estride, b, d] {
        return internals::executeInitializeFunction(m, ifct, inputs, estride, b,
                                                    b + d);
      }));
      b += d;
    }
    auto res = MultiThreadedBehaviourIntegrationResult{};
    for (auto& t : tasks) {
      const auto& ri = *(t.get());
      res.exit_status = std::min(res.exit_status, ri.exit_status);
      res.results.push_back(ri);
    }
    return res;
  }  // end of executeInitializeFunction

  BehaviourIntegrationResult integrate(MaterialDataManager& m,
                                       const BehaviourIntegrationOptions& opts,
                                       const real dt) {
    internals::allocate(m, opts);
    return internals::integrate(m, opts, dt, 0, m.n);
  }  // end of integrate

  BehaviourIntegrationResult integrate(MaterialDataManager& m,
                                       const BehaviourIntegrationOptions& opts,
                                       const real dt,
                                       const size_type b,
                                       const size_type e) {
    internals::allocate(m, opts);
    internals::checkIntegrationPointsRange(m, b, e);
    return internals::integrate(m, opts, dt, b, e);
  }  // end of integrate

  MultiThreadedBehaviourIntegrationResult integrate(
      mgis::ThreadPool& p,
      MaterialDataManager& m,
      const BehaviourIntegrationOptions& opts,
      const real dt) {
    m.setThreadSafe(true);
    internals::allocate(m, opts);
    // get number of threads
    const auto nth = p.getNumberOfThreads();
    const auto d = m.n / nth;
    const auto r = m.n % nth;
    size_type b = 0;
    std::vector<std::future<ThreadedTaskResult<BehaviourIntegrationResult>>>
        tasks;
    tasks.reserve(nth);
    for (size_type i = 0; i != r; ++i) {
      tasks.push_back(p.addTask([&m, &opts, dt, b, d] {
        return internals::integrate(m, opts, dt, b, b + d + 1);
      }));
      b += d + 1;
    }
    for (size_type i = r; i != nth; ++i) {
      tasks.push_back(p.addTask([&m, &opts, dt, b, d] {
        return internals::integrate(m, opts, dt, b, b + d);
      }));
      b += d;
    }
    auto res = MultiThreadedBehaviourIntegrationResult{};
    for (auto& t : tasks) {
      const auto& ri = *(t.get());
      res.exit_status = std::min(res.exit_status, ri.exit_status);
      res.results.push_back(ri);
    }
    return res;
  }  // end of integrate

  int integrate(mgis::ThreadPool& p,
                MaterialDataManager& m,
                const IntegrationType it,
                const real dt) {
    BehaviourIntegrationOptions opts;
    opts.integration_type = it;
    const auto r = integrate(p, m, opts, dt);
    return r.exit_status;
  }  // end of integrate

  int integrate(MaterialDataManager& m,
                const IntegrationType it,
                const real dt,
                const size_type b,
                const size_type e) {
    BehaviourIntegrationOptions opts;
    opts.integration_type = it;
    const auto r = integrate(m, opts, dt, b, e);
    return r.exit_status;
  }  // end of integrate

  static const BehaviourPostProcessing& getBehaviourPostProcessing(
      const Behaviour& b, const std::string_view n) {
    const auto p = b.postprocessings.find(n);
    if (p == b.postprocessings.end()) {
      mgis::raise(
          "getBehaviourPostProcessing: "
          "no postprocessing named '" +
          std::string{n} + "'");
    }
    return p->second;
  }  // end of getBehaviourPostProcessing

  int executePostProcessing(std::span<real> outputs,
                            BehaviourDataView& d,
                            const Behaviour& b,
                            const std::string_view n) {
    const auto& p = getBehaviourPostProcessing(b, n);
    if (outputs.size() != getArraySize(p.outputs, b.hypothesis)) {
      mgis::raise(
          "executePostProcessing: "
          "invalid size of the outputs '" +
          std::string{n} + "'");
    }
    return (p.f)(outputs.data(), &d);
  }  // end of executePostProcessing

  BehaviourIntegrationResult executePostProcessing(std::span<real> outputs,
                                                   MaterialDataManager& m,
                                                   const std::string_view n,
                                                   const size_type b,
                                                   const size_type e) {
    const auto& p = getBehaviourPostProcessing(m.b, n);
    const auto ostride = getArraySize(p.outputs, m.b.hypothesis);
    internals::checkIntegrationPointsRange(m, b, e);
    if (outputs.size() != m.n * ostride) {
      mgis::raise(
          "executePostProcessing: "
          "invalid size of the outputs '" +
          std::string{n} + "'");
    }
    return internals::executePostProcessing(outputs, m, p, ostride, b, e);
  }  // end of executePostProcessing

  BehaviourIntegrationResult executePostProcessing(std::span<real> outputs,
                                                   MaterialDataManager& m,
                                                   const std::string_view n) {
    return executePostProcessing(outputs, m, n, 0, m.n);
  }  // end of executePostProcessing

  MultiThreadedBehaviourIntegrationResult executePostProcessing(
      std::span<real> outputs,
      mgis::ThreadPool& p,
      MaterialDataManager& m,
      const std::string_view n) {
    const auto& post = getBehaviourPostProcessing(m.b, n);
    const auto ostride = getArraySize(post.outputs, m.b.hypothesis);
    if (outputs.size() != m.n * ostride) {
      mgis::raise(
          "executePostProcessing: "
          "invalid size of the outputs '" +
          std::string{n} + "'");
    }
    m.setThreadSafe(true);
    // get number of threads
    const auto nth = p.getNumberOfThreads();
    const auto d = m.n / nth;
    const auto r = m.n % nth;
    size_type b = 0;
    std::vector<std::future<ThreadedTaskResult<BehaviourIntegrationResult>>>
        tasks;
    tasks.reserve(nth);
    for (size_type i = 0; i != r; ++i) {
      tasks.push_back(p.addTask([&outputs, &m, &post, ostride, b, d] {
        return internals::executePostProcessing(outputs, m, post, ostride, b,
                                                b + d + 1);
      }));
      b += d + 1;
    }
    for (size_type i = r; i != nth; ++i) {
      tasks.push_back(p.addTask([&outputs, &m, &post, ostride, b, d] {
        return internals::executePostProcessing(outputs, m, post, ostride, b,
                                                b + d);
      }));
      b += d;
    }
    auto res = MultiThreadedBehaviourIntegrationResult{};
    for (auto& t : tasks) {
      const auto& ri = *(t.get());
      res.exit_status = std::min(res.exit_status, ri.exit_status);
      res.results.push_back(ri);
    }
    return res;
  }  // end of executePostProcessing

}  // end of namespace mgis::behaviour
