/*!
 * \file   StateView.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   02/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include "MGIS/Behaviour/State.hxx"
#include "MGIS/Behaviour/StateView.hxx"

namespace mgis {

  namespace behaviour {

    StateView make_view(State& s) { StateView v;
      auto get_ptr = [](const std::vector<real>& v) -> real* {
        if (v.empty()) {
          return nullptr;
        }
        return &v[0];
      };  // end of get_ptr
      v.gradients = get_ptr(s.gradients);
      v.thermodynamic_forces = get_ptr(s.thermodynamic_forces);
      v.material_properties = get_ptr(s.material_properties);
      v.internal_state_variables = get_ptr(s.internal_state_variables);
      v.external_state_variables = get_ptr(s.external_state_variables);
      return v;
    }  // end of make_view

  }  // end of namepace behaviour

}  // end of namespace mgis