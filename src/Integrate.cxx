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

#include <tuple>
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

namespace mgis {

  namespace behaviour {

    int integrate(MaterialDataManager& m,
                  const size_type b,
                  const size_type e) {
      /*
       * \brief uniform values are treated immediatly. For spatially variable
       * fields, we return the information needed to evaluate them
       */
      auto dispatch = [](
          std::vector<real>& v,
          std::map<std::string, mgis::variant<real, std::vector<real>>>&
              values) {
        std::vector<std::tuple<size_type, real*>> evs;
        size_type i = 0;
        for (auto& value : values) {
          if (holds_alternative<real>(value.second)) {
            v[i] = get<real>(value.second);
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
      const auto g_stride = getArraySize(m.b.gradients, m.b.hypothesis);
      const auto t_stride =
          getArraySize(m.b.thermodynamic_forces, m.b.hypothesis);
      const auto isvs_stride = getArraySize(m.b.isvs, m.b.hypothesis);
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
        eval(mps0, vmps0, i);
        eval(mps1, vmps1, i);
        eval(esvs0, vesvs0, i);
        eval(esvs1, vesvs1, i);
        v.s0.gradients = m.s0.gradients.data() + g_stride;
        v.s1.gradients = m.s1.gradients.data() + g_stride;
        v.s0.thermodynamic_forces =
            m.s0.thermodynamic_forces.data() + t_stride;
        v.s1.thermodynamic_forces =
            m.s1.thermodynamic_forces.data() + t_stride;
        v.s0.material_properties = mps0.data();
        v.s1.material_properties = mps1.data();
        v.s0.internal_state_variables =
            m.s0.internal_state_variables.data() + isvs_stride;
        v.s1.internal_state_variables =
            m.s1.internal_state_variables.data() + isvs_stride;
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

  }  // end of namespace behaviour

}  // end of namespace mgis
