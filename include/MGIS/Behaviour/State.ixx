/*!
 * \file   State.ixx
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

#ifndef LIB_MGIS_BEHAVIOUR_STATE_IXX
#define LIB_MGIS_BEHAVIOUR_STATE_IXX

namespace mgis::behaviour {

  inline StateView make_view(State& s) {
    auto get_ptr = [](std::vector<real>& v) -> real* {
      if (v.empty()) {
        return nullptr;
      }
      return &v[0];
    };  // end of get_ptr
    StateView v;
    v.mass_density = &(s.mass_density);
    v.gradients = get_ptr(s.gradients);
    v.thermodynamic_forces = get_ptr(s.thermodynamic_forces);
    v.material_properties = get_ptr(s.material_properties);
    v.internal_state_variables = get_ptr(s.internal_state_variables);
    v.stored_energy = &s.stored_energy;
    v.dissipated_energy = &s.dissipated_energy;
    v.external_state_variables = get_ptr(s.external_state_variables);
    return v;
  }  // end of make_view

  inline InitialStateView make_view(const State& s) {
    auto get_ptr = [](const std::vector<real>& v) -> const real* {
      if (v.empty()) {
        return nullptr;
      }
      return &v[0];
    };  // end of get_ptr
    InitialStateView v;
    v.mass_density = &(s.mass_density);
    v.gradients = get_ptr(s.gradients);
    v.thermodynamic_forces = get_ptr(s.thermodynamic_forces);
    v.material_properties = get_ptr(s.material_properties);
    v.internal_state_variables = get_ptr(s.internal_state_variables);
    v.stored_energy = &s.stored_energy;
    v.dissipated_energy = &s.dissipated_energy;
    v.external_state_variables = get_ptr(s.external_state_variables);
    return v;
  }  // end of make_view

}  // end of namespace mgis::behaviour

#endif /* LIB_MGIS_BEHAVIOUR_STATE_IXX */
