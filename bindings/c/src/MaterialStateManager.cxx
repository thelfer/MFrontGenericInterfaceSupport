/*!
 * \file   bindings/c/src/MaterialStateManager.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   11/09/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include "MGIS/Behaviour/MaterialStateManager.h"

extern "C" {

mgis_status mgis_bv_material_state_manager_get_gradients(
    mgis_real** g, mgis_bv_MaterialStateManager* const s) {
  if (s == nullptr) {
    *g = nullptr;
    return mgis_report_failure("null state manager");
  }
  *g = s->gradients.data();
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_get_gradients

mgis_status mgis_bv_material_state_manager_get_gradients_stride(
    mgis_size_type* const gs, mgis_bv_MaterialStateManager* const s) {
  if (s == nullptr) {
    *gs = mgis_size_type{};
    return mgis_report_failure("null state manager");
  }
  *gs = s->gradients_stride;
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_get_gradients_stride

mgis_status mgis_bv_material_state_manager_get_thermodynamic_forces(
    mgis_real** g, mgis_bv_MaterialStateManager* const s) {
  if (s == nullptr) {
    *g = nullptr;
    return mgis_report_failure("null state manager");
  }
  *g = s->thermodynamic_forces.data();
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_get_thermodynamic_forces

mgis_status mgis_bv_material_state_manager_get_thermodynamic_forces_stride(
    mgis_size_type* const gs, mgis_bv_MaterialStateManager* const s) {
  if (s == nullptr) {
    *gs = mgis_size_type{};
    return mgis_report_failure("null state manager");
  }
  *gs = s->thermodynamic_forces_stride;
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_get_thermodynamic_forces_stride

mgis_status mgis_bv_material_state_manager_get_internal_state_variables(
    mgis_real** g, mgis_bv_MaterialStateManager* const s) {
  if (s == nullptr) {
    *g = nullptr;
    return mgis_report_failure("null state manager");
  }
  *g = s->internal_state_variables.data();
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_get_internal_state_variables

mgis_status mgis_bv_material_state_manager_get_internal_state_variables_stride(
    mgis_size_type* const gs, mgis_bv_MaterialStateManager* const s) {
  if (s == nullptr) {
    *gs = mgis_size_type{};
    return mgis_report_failure("null state manager");
  }
  *gs = s->internal_state_variables_stride;
  return mgis_report_success();
}  // end of mgis_bv_material_state_manager_get_internal_state_variables_stride


}  // end of extern "C"