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

namespace mgis {

  namespace behaviour {

    IntegrateWorkSpace::IntegrateWorkSpace(const Behaviour& b)
        : mps0(getArraySize(b.mps, b.hypothesis)),
          mps1(getArraySize(b.mps, b.hypothesis)),
          esvs0(getArraySize(b.esvs, b.hypothesis)),
          esvs1(getArraySize(b.esvs, b.hypothesis)) {
    }  // end of IntegrateWorkSpace

    IntegrateWorkSpace::IntegrateWorkSpace(IntegrateWorkSpace&&) = default;
    IntegrateWorkSpace::IntegrateWorkSpace(const IntegrateWorkSpace&) = default;
    IntegrateWorkSpace& IntegrateWorkSpace::operator=(IntegrateWorkSpace&&) =
        default;
    IntegrateWorkSpace& IntegrateWorkSpace::operator=(
        const IntegrateWorkSpace&) = default;

    IntegrateWorkSpace& getIntegrateWorkSpace(const Behaviour& b) {
      static std::map<
          const Behaviour*,
          std::map<std::thread::id, std::shared_ptr<IntegrateWorkSpace>>>
          m;
      static std::mutex mt;
      const auto id = std::this_thread::get_id();
      std::lock_guard<std::mutex> lock(mt);
      auto& mwks = m[&b];
      auto p = mwks.find(id);
      if (p == mwks.end()) {
        p = mwks.insert({id, std::make_shared<IntegrateWorkSpace>(b)}).first;
      }
      return *(p->second);
    }  // end of getIntegrateWorkSpace

    int integrate(MaterialDataManager& m,
                  const IntegrationType it,
                  const real dt,
                  const size_type b,
                  const size_type e) {
      /*
       * \brief uniform values are treated immediatly. For spatially variable
       * fields, we return the information needed to evaluate them
       */
      auto dispatch = [](
          std::vector<real>& v,
          std::map<std::string, mgis::variant<real, mgis::span<real>,
                                              std::vector<real>>>& values,
          const std::vector<Variable>& ds) {
        mgis::raise_if(ds.size() != v.size(),
                       "integrate: ill allocated memory");
        // evaluators
        std::vector<std::tuple<size_type, real*>> evs;
        auto i = mgis::size_type{};
        for (const auto& d : ds) {
          if (d.type != Variable::SCALAR) {
            mgis::raise("integrate: invalid type for variable '" + d.name +
                        "'");
          }
          auto p = values.find(d.name);
          if (p == values.end()) {
            auto msg = std::string{"integrate: no variable named '" + d.name +
                                   "' declared"};
            if (!values.empty()) {
              msg += "\nThe following variables were declared: ";
              for (const auto& vs : values) {
                msg += "\n- " + vs.first;
              }
            } else {
              msg += "\nNo variable declared.";
            }
            mgis::raise(msg);
          }
          if (holds_alternative<real>(p->second)) {
            v[i] = get<real>(p->second);
          } else if (holds_alternative<mgis::span<real>>(p->second)) {
            evs.push_back(
                std::make_tuple(i, get<mgis::span<real>>(p->second).data()));
          } else {
            evs.push_back(
                std::make_tuple(i, get<std::vector<real>>(p->second).data()));
          }
          ++i;
        }
        return evs;
      };  // end of dispatch
      auto eval = [](std::vector<real>& v,
                     const std::vector<std::tuple<size_type, real*>>& evs,
                     const size_type i) {
        for (const auto& ev : evs) {
          v[std::get<0>(ev)] = std::get<1>(ev)[i];
        }
      };  // end of eval
      // strides
      const auto g_stride = m.s0.gradients_stride;
      const auto t_stride = m.s0.thermodynamic_forces_stride;
      const auto isvs_stride = m.s0.internal_state_variables_stride;
      // workspace
      auto& ws = getIntegrateWorkSpace(m.b);
      // treating uniform values
      const auto vmps0 = dispatch(ws.mps0, m.s0.material_properties, m.b.mps);
      const auto vmps1 = dispatch(ws.mps1, m.s1.material_properties, m.b.mps);
      const auto vesvs0 =
          dispatch(ws.esvs0, m.s0.external_state_variables, m.b.esvs);
      const auto vesvs1 =
          dispatch(ws.esvs1, m.s1.external_state_variables, m.b.esvs);
      // loop over integration points
      auto r = int{1};
      real opts[Behaviour::nopts + 1];  // option passed to the behaviour
      for (auto i = b; i != e; ++i) {
        BehaviourDataView v;
        v.rdt = real(1);
        v.dt = dt;
        if (it != IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR) {
          v.K = m.K.data() + m.K_stride * i;
        } else {
          v.K = &opts[0];
        }
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
        v.s0.stored_energy = m.s0.stored_energies.data() + i;
        v.s1.stored_energy = m.s1.stored_energies.data() + i;
        v.s0.dissipated_energy = m.s0.dissipated_energies.data() + i;
        v.s1.dissipated_energy = m.s1.dissipated_energies.data() + i;
        v.s0.external_state_variables = ws.esvs0.data();
        v.s1.external_state_variables = ws.esvs1.data();
        v.K[0] = static_cast<int>(it);
        switch (integrate(v, m.b)) {
          case 1:
            r = std::min(r, 1);
            break;
          case 0:
            r = 0;
            break;
          default:
            return -1;
        }
      }
      return r;
    }  // end of integrate

    int integrate(ThreadPool& p,
                  MaterialDataManager& m,
                  const IntegrationType it,
                  const real dt) {
      // get number of threads
      const auto nth = p.getNumberOfThreads();
      const auto d = m.n / nth;
      const auto r = m.n % nth;
      size_type b = 0;
      std::vector<std::future<ThreadedTaskResult<int>>> tasks;
      tasks.reserve(nth);
      for (size_type i = 0; i != r; ++i) {
        tasks.push_back(p.addTask(
            [&m, it, dt, b, d] { return integrate(m, it, dt, b, b + d + 1); }));
        b += d + 1;
      }
      for (size_type i = r; i != nth; ++i) {
        tasks.push_back(p.addTask(
            [&m, it, dt, b, d] { return integrate(m, it, dt, b, b + d); }));
        b += d;
      }
      //      p.wait();
      auto res = int{1};
      for (auto& t : tasks) {
        switch (*(t.get())) {
          case 1:
            res = std::min(res, 1);
            break;
          case 0:
            res = 0;
            break;
          default:
            return -1;
        }
      }
      return res;
    }  // end of integrate

  }  // end of namespace behaviour

}  // end of namespace mgis
