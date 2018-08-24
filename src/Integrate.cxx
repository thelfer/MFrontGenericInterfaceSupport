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

#include<iostream>

#include <tuple>
#include <cstdlib>
#include <cinttypes>
#include "MGIS/ThreadPool.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

namespace mgis {

  namespace behaviour {

    int integrate(MaterialDataManager& m,
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
                                              std::vector<real>>>& values) {
        std::vector<std::tuple<size_type, real*>> evs;
        size_type i = 0;
        for (auto& value : values) {
          if (holds_alternative<real>(value.second)) {
            v[i] = get<real>(value.second);
          } else if (holds_alternative<mgis::span<real>>(value.second)) {
            evs.push_back(std::make_tuple(
                i, get<mgis::span<real>>(value.second).data()));
          } else {
            evs.push_back(std::make_tuple(
                i, get<std::vector<real>>(value.second).data()));
          }
          ++i;
        }
        return evs;
      };  // end of dispatch
      auto eval = [](std::vector<real>& v,
                     const std::vector<std::tuple<size_type, real*>>& evs,
                     const size_type i) {
        for (const auto ev : evs) {
          v[std::get<0>(ev)] = std::get<1>(ev)[i];
        }
      };  // end of eval
      // strides
      const auto g_stride = m.s0.gradients_stride;
      const auto t_stride = m.s0.thermodynamic_forces_stride;
      const auto isvs_stride = m.s0.internal_state_variables_stride;
      // workspace
      std::vector<real> mps0(getArraySize(m.b.mps, m.b.hypothesis));
      std::vector<real> mps1(getArraySize(m.b.mps, m.b.hypothesis));
      std::vector<real> esvs0(getArraySize(m.b.esvs, m.b.hypothesis));
      std::vector<real> esvs1(getArraySize(m.b.esvs, m.b.hypothesis));
      // treating uniform values
      const auto vmps0 = dispatch(mps0, m.s0.material_properties);
      const auto vmps1 = dispatch(mps1, m.s1.material_properties);
      const auto vesvs0 = dispatch(esvs0, m.s0.external_state_variables);
      const auto vesvs1 = dispatch(esvs1, m.s1.external_state_variables);
      // loop over integration points
      auto r = int{1};
      for (auto i = b; i != e + 1; ++i) {
        BehaviourDataView v;
        v.dt = dt;
        v.K = m.K.data() + m.K_stride * i;
        eval(mps0, vmps0, i);
        eval(mps1, vmps1, i);
        eval(esvs0, vesvs0, i);
        eval(esvs1, vesvs1, i);
        v.s0.gradients = m.s0.gradients.data() + g_stride * i;
        v.s1.gradients = m.s1.gradients.data() + g_stride * i;
        v.s0.thermodynamic_forces =
            m.s0.thermodynamic_forces.data() + t_stride * i;
        v.s1.thermodynamic_forces =
            m.s1.thermodynamic_forces.data() + t_stride * i;
        v.s0.material_properties = mps0.data();
        v.s1.material_properties = mps1.data();
        v.s0.internal_state_variables =
            m.s0.internal_state_variables.data() + isvs_stride * i;
        v.s1.internal_state_variables =
            m.s1.internal_state_variables.data() + isvs_stride * i;
        v.s0.external_state_variables = esvs0.data();
        v.s1.external_state_variables = esvs1.data();
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

    int integrate(ThreadPool& p, MaterialDataManager& m, const real dt) {
      // get number of threads
      const auto nth = p.getNumberOfThreads();
      const auto d = m.n / nth;
      const auto r = m.n % nth;
      size_type b = 0;
      std::vector<std::future<ThreadedTaskResult<int>>> tasks;
      tasks.reserve(nth);
      for (size_type i = 0; i != r; ++i) {
        tasks.push_back(
            p.addTask([&m, dt, b, d] { return integrate(m, dt, b, b+d+1); }));
        b += d+1;
      }
      for (size_type i = r; i != nth; ++i) {
        tasks.push_back(
            p.addTask([&m, dt, b, d] { return integrate(m, dt, b, b+d); }));
        b += d ;
      }
      p.wait();
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
