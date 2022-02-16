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

namespace mgis::behaviour {

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

  /*!
   * \brief integration the behaviour over a range of integration points.
   */
  static BehaviourIntegrationResult integrate2(
      MaterialDataManager& m,
      const BehaviourIntegrationOptions& opts,
      const real dt,
      const size_type b,
      const size_type e) {
    /*
     * \brief uniform values are treated immediatly. For spatially variable
     * fields, we return the information needed to evaluate them
     */
    auto dispatch =
        [&m](std::vector<real>& v,
             std::map<std::string, MaterialStateManager::FieldHolder>& values,
             const std::vector<Variable>& ds) {
          mgis::raise_if(v.size() != getArraySize(ds, m.b.hypothesis),
                         "integrate: ill allocated memory");
          // evaluators
          std::vector<std::tuple<size_type,  // offset
                                 size_type,  // variable size
                                 const real*>>
              evs;
          auto offset = mgis::size_type{};
          for (const auto& d : ds) {
            const auto s = getVariableSize(d, m.b.hypothesis);
            auto p = values.find(d.name);
            if (p == values.end()) {
              auto msg = std::string{"integrate: no variable named '" + d.name +
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
                evs.push_back(
                    std::make_tuple(offset, s, variable_values.data()));
              }
            };
            if (std::holds_alternative<real>(p->second)) {
              if (d.type != Variable::SCALAR) {
                mgis::raise("integrate: invalid type for variable '" + d.name +
                            "'");
              }
              v[offset] = std::get<real>(p->second);
            } else if (std::holds_alternative<mgis::span<real>>(p->second)) {
              set_values(std::get<mgis::span<real>>(p->second));
            } else {
              set_values(std::get<std::vector<real>>(p->second));
            }
            offset += s;
          }
          return evs;
        };  // end of dispatch
    auto eval =
        [](std::vector<real>& values,
           const std::vector<std::tuple<size_type, size_type, const real*>>& evs,
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
        };  // end of eval
    // strides
    const auto g_stride = m.s0.gradients_stride;
    const auto t_stride = m.s0.thermodynamic_forces_stride;
    const auto isvs_stride = m.s0.internal_state_variables_stride;
    // workspace
    auto& ws = m.getBehaviourIntegrationWorkSpace();
    // treating uniform values
    const auto vmps0 = dispatch(ws.mps0, m.s0.material_properties, m.b.mps);
    const auto vmps1 = dispatch(ws.mps1, m.s1.material_properties, m.b.mps);
    const auto vesvs0 =
        dispatch(ws.esvs0, m.s0.external_state_variables, m.b.esvs);
    const auto vesvs1 =
        dispatch(ws.esvs1, m.s1.external_state_variables, m.b.esvs);
    const auto computes_stored_energy = m.b.computesStoredEnergy;
    const auto computes_dissipated_energy = m.b.computesDissipatedEnergy;
    // loop over integration points
    auto r = BehaviourIntegrationResult{};
    auto rdt0 = r.time_step_increase_factor;
    const real Ke = encodeBehaviourIntegrationOptions(opts);
    real bopts[Behaviour::nopts + 1];  // option passed to the behaviour
    for (auto i = b; i != e; ++i) {
      auto rdt = rdt0;
      BehaviourDataView v;
      v.error_message = ws.error_message.data();
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
      v.speed_of_sound = m.speed_of_sound.data() + i;
      eval(ws.mps0, vmps0, i);
      eval(ws.mps1, vmps1, i);
      eval(ws.esvs0, vesvs0, i);
      eval(ws.esvs1, vesvs1, i);
      v.s0.gradients = m.s0.gradients.data() + g_stride * i;
      v.s1.gradients = m.s1.gradients.data() + g_stride * i;
      v.s0.thermodynamic_forces =
          m.s0.thermodynamic_forces.data() + t_stride * i;
      v.s1.thermodynamic_forces =
          m.s1.thermodynamic_forces.data() + t_stride * i;
      v.s0.material_properties = ws.mps0.data();
      v.s1.material_properties = ws.mps1.data();
      v.s0.internal_state_variables =
          m.s0.internal_state_variables.data() + isvs_stride * i;
      v.s1.internal_state_variables =
          m.s1.internal_state_variables.data() + isvs_stride * i;
      if (computes_stored_energy) {
        v.s0.stored_energy = m.s0.stored_energies.data() + i;
        v.s1.stored_energy = m.s1.stored_energies.data() + i;
      } else {
        v.s0.stored_energy = nullptr;
        v.s1.stored_energy = nullptr;
      }
      if (computes_dissipated_energy) {
        v.s0.dissipated_energy = m.s0.dissipated_energies.data() + i;
        v.s1.dissipated_energy = m.s1.dissipated_energies.data() + i;
      } else {
        v.s0.dissipated_energy = nullptr;
        v.s1.dissipated_energy = nullptr;
      }
      v.s0.external_state_variables = ws.esvs0.data();
      v.s1.external_state_variables = ws.esvs1.data();
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
  }  // end of integrate2

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

  BehaviourIntegrationResult integrate(MaterialDataManager& m,
                                       const BehaviourIntegrationOptions& opts,
                                       const real dt,
                                       const size_type b,
                                       const size_type e) {
    allocate(m, opts);
    return integrate2(m, opts, dt, b, e);
  }  // end of integrate

  int integrate(ThreadPool& p,
                MaterialDataManager& m,
                const IntegrationType it,
                const real dt) {
    BehaviourIntegrationOptions opts;
    opts.integration_type = it;
    const auto r = integrate(p, m, opts, dt);
    return r.exit_status;
  }  // end of integrate

  MultiThreadedBehaviourIntegrationResult integrate(
      ThreadPool& p,
      MaterialDataManager& m,
      const BehaviourIntegrationOptions& opts,
      const real dt) {
    m.setThreadSafe(true);
    allocate(m, opts);
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
        return integrate2(m, opts, dt, b, b + d + 1);
      }));
      b += d + 1;
    }
    for (size_type i = r; i != nth; ++i) {
      tasks.push_back(p.addTask(
          [&m, &opts, dt, b, d] { return integrate2(m, opts, dt, b, b + d); }));
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

  int executeInitializeFunction(BehaviourDataView& d,
                                const Behaviour& b,
                                const std::string_view n) {
    const auto p = b.initialize_functions.find(n);
    if (p == b.initialize_functions.end()) {
      mgis::raise(
          "executeInitializeFunction: "
          "no initialize function named '" +
          std::string{n} + "'");
    }
    if (!p->second.inputs.empty()) {
      mgis::raise(
          "executeInitializeFunction: "
          "invalid size of the inputs '" +
          std::string{n} + "'");
    }
    const auto& f = p->second.f;
    return f(&d, nullptr);
  }  // end of executeInitializeFunction


  int executeInitializeFunction(BehaviourDataView& d,
                                mgis::span<const real> inputs,
                                const Behaviour& b,
                                const std::string_view n) {
    const auto p = b.initialize_functions.find(n);
    if (p == b.initialize_functions.end()) {
      mgis::raise(
          "executeInitializeFunction: "
          "no initialize function named '" +
          std::string{n} + "'");
    }
    if (inputs.size() != getArraySize(p->second.inputs, b.hypothesis)) {
      mgis::raise(
          "executeInitializeFunction: "
          "invalid size of the inputs '" +
          std::string{n} + "'");
    }
    const auto& f = p->second.f;
    return f(&d, inputs.data());
  }  // end of executeInitializeFunction

  int executePostProcessing(mgis::span<real> outputs,
                            BehaviourDataView& d,
                            const Behaviour& b,
                            const std::string_view n) {
    const auto p = b.postprocessings.find(n);
    if (p == b.postprocessings.end()) {
      mgis::raise(
          "executePostProcessing: "
          "no postprocessing named '" +
          std::string{n} + "'");
    }
    if (outputs.size() != getArraySize(p->second.outputs, b.hypothesis)) {
      mgis::raise(
          "executePostProcessing: "
          "invalid size of the outputs '" +
          std::string{n} + "'");
    }
    const auto& f = p->second.f;
    return f(outputs.data(), &d);
  }  // end of executePostProcessing

}  // end of namespace mgis::behaviour
