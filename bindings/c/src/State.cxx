/*!
 * \file   bindings/c/src/State.cxx
 * \brief
 * \author Thomas Helfer
 * \date   03/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include "MGIS/Behaviour/Variable.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/State.h"

extern "C" {

mgis_status mgis_bv_state_set_mass_density(mgis_bv_State* const s,
                                           const mgis_real v) {
  if (s == nullptr) {
    return mgis_report_failure("invalid argument (null state)");
  }
  s->mass_density = v;
  return mgis_report_success();
}  // end of mgis_bv_state_set_mass_density

mgis_status mgis_bv_state_set_gradient_by_name(mgis_bv_State* const s,
                                               const char* const n,
                                               const mgis_real* const v) {
  if (s == nullptr) {
    return mgis_report_failure("invalid argument (null state)");
  }
  if (n == nullptr) {
    return mgis_report_failure("invalid argument (null name)");
  }
  if (v == nullptr) {
    return mgis_report_failure("invalid argument (null values)");
  }
  try {
    setGradient(*s, n, v);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_set_gradient_by_name

mgis_status mgis_bv_state_get_gradient_by_name(mgis_real** v,
                                               mgis_bv_State* const s,
                                               const char* const n) {
  if (s == nullptr) {
    return mgis_report_failure("invalid argument (null state)");
  }
  if (n == nullptr) {
    return mgis_report_failure("invalid argument (null name)");
  }
  if (v == nullptr) {
    return mgis_report_failure("invalid argument (null values)");
  }
  try {
    *v = getGradient(*s, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_get_gradient_by_name

mgis_status mgis_bv_state_set_gradient_by_offset(mgis_bv_State* const s,
                                                 const mgis_size_type o,
                                                 const mgis_size_type n,
                                                 const mgis_real* const v) {
  if (s == nullptr) {
    return mgis_report_failure("invalid argument (null state)");
  }
  if (v == nullptr) {
    return mgis_report_failure("invalid argument (null values)");
  }
  if (n == 0) {
    return mgis_report_failure("invalid argument (null size)");
  }
  try {
    setGradient(*s, o, n, v);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_set_gradient_by_offset

mgis_status mgis_bv_state_get_gradient_by_offset(mgis_real** v,
                                                 mgis_bv_State* const s,
                                                 const mgis_size_type o) {
  if (v == nullptr) {
    return mgis_report_failure("invalid argument (null values)");
  }
  if (s == nullptr) {
    return mgis_report_failure("invalid argument (null state)");
  }
  try {
    *v = getGradient(*s, o);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_get_gradient_by_offset

mgis_status mgis_bv_state_set_thermodynamic_force_by_name(
    mgis_bv_State* const s, const char* const n, const mgis_real* const v) {
  if (s == nullptr) {
    return mgis_report_failure("invalid argument (null state)");
  }
  if (n == nullptr) {
    return mgis_report_failure("invalid argument (null name)");
  }
  if (v == nullptr) {
    return mgis_report_failure("invalid argument (null values)");
  }
  try {
    setThermodynamicForce(*s, n, v);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_set_thermodynamic_force_by_name

mgis_status mgis_bv_state_get_thermodynamic_force_by_name(
    mgis_real** v, mgis_bv_State* const s, const char* const n) {
  if (s == nullptr) {
    return mgis_report_failure("invalid argument (null state)");
  }
  if (n == nullptr) {
    return mgis_report_failure("invalid argument (null name)");
  }
  if (v == nullptr) {
    return mgis_report_failure("invalid argument (null values)");
  }
  try {
    *v = getThermodynamicForce(*s, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_get_thermodynamic_force_by_name

mgis_status mgis_bv_state_set_thermodynamic_force_by_offset(
    mgis_bv_State* const s,
    const mgis_size_type o,
    const mgis_size_type n,
    const mgis_real* const v) {
  if (s == nullptr) {
    return mgis_report_failure("invalid argument (null state)");
  }
  if (v == nullptr) {
    return mgis_report_failure("invalid argument (null values)");
  }
  if (n == 0) {
    return mgis_report_failure("invalid argument (null size)");
  }
  try {
    setThermodynamicForce(*s, o, n, v);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_set_thermodynamic_force_by_offset

mgis_status mgis_bv_state_get_thermodynamic_force_by_offset(
    mgis_real** v, mgis_bv_State* const s, const mgis_size_type o) {
  if (v == nullptr) {
    return mgis_report_failure("invalid argument (null values)");
  }
  if (s == nullptr) {
    return mgis_report_failure("invalid argument (null state)");
  }
  try {
    *v = getThermodynamicForce(*s, o);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_get_thermodynamic_force_by_offset

mgis_status mgis_bv_state_set_material_property_by_name(mgis_bv_State* const s,
                                                        const char* const n,
                                                        const mgis_real v) {
  try {
    setMaterialProperty(*s, n, v);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_set_material_property_by_name

mgis_status mgis_bv_state_get_material_property_by_name(mgis_real** const v,
                                                        mgis_bv_State* const s,
                                                        const char* const n) {
  try {
    *v = getMaterialProperty(*s, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_get_material_property_by_name

mgis_status mgis_bv_state_set_material_property_by_offset(
    mgis_bv_State* const s, const mgis_size_type o, const mgis_real v) {
  try {
    setMaterialProperty(*s, o, v);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_set_material_property_by_offset

mgis_status mgis_bv_state_get_material_property_by_offset(
    mgis_real** const v, mgis_bv_State* const s, const mgis_size_type o) {
  try {
    *v = getMaterialProperty(*s, o);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_get_material_property_by_offset
  
mgis_status mgis_bv_state_set_internal_state_variable_by_name(
    mgis_bv_State* const s, const char* const n, const mgis_real* const v) {
  if (s == nullptr) {
    return mgis_report_failure("invalid argument (null state)");
  }
  if (n == nullptr) {
    return mgis_report_failure("invalid argument (null name)");
  }
  if (v == nullptr) {
    return mgis_report_failure("invalid argument (null values)");
  }
  try {
    setInternalStateVariable(*s, n, v);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_set_internal_state_variable_by_name

mgis_status mgis_bv_state_get_internal_state_variables(
    mgis_real** v, mgis_bv_State* const s) {
  if (s == nullptr) {
    return mgis_report_failure("invalid argument (null state)");
  }
  if (v == nullptr) {
    return mgis_report_failure("invalid argument (null values)");
  }
  if(s->internal_state_variables.empty()){
    *v = nullptr;
    return mgis_report_failure("no internal state variables declared");
  }
  *v = s->internal_state_variables.data();
  return mgis_report_success();
}  // end of mgis_bv_state_get_internal_state_variable_by_name
  
mgis_status mgis_bv_state_get_internal_state_variable_by_name(
    mgis_real** v, mgis_bv_State* const s, const char* const n) {
  if (s == nullptr) {
    return mgis_report_failure("invalid argument (null state)");
  }
  if (n == nullptr) {
    return mgis_report_failure("invalid argument (null name)");
  }
  if (v == nullptr) {
    return mgis_report_failure("invalid argument (null values)");
  }
  try {
    *v = getInternalStateVariable(*s, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_get_internal_state_variable_by_name

mgis_status mgis_bv_state_set_internal_state_variable_by_offset(
    mgis_bv_State* const s,
    const mgis_size_type o,
    const mgis_size_type n,
    const mgis_real* const v) {
  if (s == nullptr) {
    return mgis_report_failure("invalid argument (null state)");
  }
  if (v == nullptr) {
    return mgis_report_failure("invalid argument (null values)");
  }
  if (n == 0) {
    return mgis_report_failure("invalid argument (null size)");
  }
  try {
    setInternalStateVariable(*s, o, n, v);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_set_internal_state_variable_by_offset

mgis_status mgis_bv_state_get_internal_state_variable_by_offset(
    mgis_real** v, mgis_bv_State* const s, const mgis_size_type o) {
  if (v == nullptr) {
    return mgis_report_failure("invalid argument (null values)");
  }
  if (s == nullptr) {
    return mgis_report_failure("invalid argument (null state)");
  }
  try {
    *v = getInternalStateVariable(*s, o);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_get_internal_state_variable_by_offset

mgis_status mgis_bv_state_set_scalar_external_state_variable_by_name(
    mgis_bv_State* const s, const char* const n, const mgis_real v) {
  try {
    setExternalStateVariable(*s, n, v);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_set_external_state_variable_by_name

mgis_status mgis_bv_state_set_external_state_variable_by_name(
    mgis_bv_State* const s, const char* const n, const mgis_real* const v) {
  try {
    const auto& ev = getVariable(s->b.esvs, n);
    const auto es = getVariableSize(ev, s->b.hypothesis);
    setExternalStateVariable(*s, n, mgis::span<const mgis::real>(v, es));
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_set_external_state_variable_by_name

mgis_status mgis_bv_state_get_external_state_variables(
    mgis_real** v, mgis_bv_State* const s) {
  if (s == nullptr) {
    return mgis_report_failure("invalid argument (null state)");
  }
  if (v == nullptr) {
    return mgis_report_failure("invalid argument (null values)");
  }
  if(s->external_state_variables.empty()){
    *v = nullptr;
    return mgis_report_failure("no external state variables declared");
  }
  *v = s->external_state_variables.data();
  return mgis_report_success();
}  // end of mgis_bv_state_get_external_state_variable_by_name

mgis_status mgis_bv_state_get_external_state_variable_by_name(
    mgis_real** const v, mgis_bv_State* const s, const char* const n) {
  try {
    *v = getExternalStateVariable(*s, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_get_external_state_variable_by_name

mgis_status mgis_bv_state_set_scalar_external_state_variable_by_offset(
    mgis_bv_State* const s, const mgis_size_type o, const mgis_real v) {
  try {
    setExternalStateVariable(*s, o, v);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_set_external_state_variable_by_offset

mgis_status mgis_bv_state_set_external_state_variable_by_offset(
    mgis_bv_State* const s,
    const mgis_size_type o,
    const mgis_real* const v,
    const mgis_size_type vs) {
  try {
    setExternalStateVariable(*s, o, mgis::span<const mgis::real>(v, vs));
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_set_external_state_variable_by_offset

mgis_status mgis_bv_state_get_external_state_variable_by_offset(
    mgis_real** const v, mgis_bv_State* const s, const mgis_size_type o) {
  try {
    *v = getExternalStateVariable(*s, o);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_state_get_external_state_variable_by_offset

}  // end of extern "C"
