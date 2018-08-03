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

#include "MGIS/Behaviour/State.h"

extern "C" {

mgis_status mgis_bv_set_state_material_property_by_name(
    mgis_bv_State* const s,
    const mgis_bv_Behaviour* const b,
    const char* const n,
    const mgis_real v) {
  try {
    setMaterialProperty(*s, *b, n, v);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_set_state_material_property_by_name

mgis_status mgis_bv_get_state_material_property_by_name(
    mgis_real* const v,
    const mgis_bv_State* const s,
    const mgis_bv_Behaviour* const b,
    const char* const n) {
  try {
    *v = getMaterialProperty(*s, *b, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_get_state_material_property_by_name

mgis_status mgis_bv_set_state_material_property_by_offset(
    mgis_bv_State* const s, const mgis_size_type o, const mgis_real v) {
  try {
    setMaterialProperty(*s, o, v);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_set_state_material_property_by_offset

mgis_status mgis_bv_get_state_material_property_by_offset(
    mgis_real* const v, const mgis_bv_State* const s, const mgis_size_type o) {
  try {
    *v = getMaterialProperty(*s, o);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_get_state_material_property_by_offset

mgis_status mgis_bv_set_state_external_state_variable_by_name(
    mgis_bv_State* const s,
    const mgis_bv_Behaviour* const b,
    const char* const n,
    const mgis_real v) {
  try {
    setExternalStateVariable(*s, *b, n, v);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_set_state_external_state_variable_by_name

mgis_status mgis_bv_get_state_external_state_variable_by_name(
    mgis_real* const v,
    const mgis_bv_State* const s,
    const mgis_bv_Behaviour* const b,
    const char* const n) {
  try {
    *v = getExternalStateVariable(*s, *b, n);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_get_state_external_state_variable_by_name

mgis_status mgis_bv_set_state_external_state_variable_by_offset(
    mgis_bv_State* const s, const mgis_size_type o, const mgis_real v) {
  try {
    setExternalStateVariable(*s,o, v);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_set_state_external_state_variable_by_offset

mgis_status mgis_bv_get_state_external_state_variable_by_offset(
    mgis_real* const v, const mgis_bv_State* const s, const mgis_size_type o) {
  try {
    *v = getExternalStateVariable(*s, o);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_get_state_external_state_variable_by_offset

} // end of extern "C"