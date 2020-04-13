/*!
 * \file   MaterialDataManager.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   05/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include "MGIS/Behaviour/MaterialDataManager.h"

extern "C" {

mgis_status mgis_bv_create_material_data_manager_initializer(
    mgis_bv_MaterialDataManagerInitializer** d) {
  *d = nullptr;
  try {
    *d = new mgis::behaviour::MaterialDataManagerInitializer();
    if (*d == nullptr) {
      return mgis_report_failure(
          "mgis_bv_create_material_data_manager_initializer: "
          "memory allocation failed");
    }
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_create_material_data_manager_initializer

mgis_status mgis_bv_material_data_manager_initializer_bind_tangent_operator(
    mgis_bv_MaterialDataManagerInitializer* d,
    mgis_real* const K,
    mgis_size_type s) {
  if (d == nullptr) {
    return mgis_report_failure(
        "mgis_bv_material_data_manager_initializer_bind_tangent_operator: "
        "null argument");
  }
  if (K == nullptr) {
    return mgis_report_failure(
        "mgis_bv_material_data_manager_initializer_bind_tangent_operator: "
        "invalid tangent operator");
  }
  try {
    d->K = mgis::span<mgis::real>(K, s);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_create_material_data_manager_initializer

mgis_status mgis_bv_material_data_manager_initializer_get_state_0_initializer(
    mgis_bv_MaterialStateManagerInitializer** s,
    mgis_bv_MaterialDataManagerInitializer* const d) {
  if (d == nullptr) {
    return mgis_report_failure("invalid argument (material data manager is null)");
  }
  *s = &(d->s0);
  return mgis_report_success();
}  // end of mgis_bv_material_data_manager_get_state_0_initializer

mgis_status mgis_bv_material_data_manager_initializer_get_state_1_initializer(
    mgis_bv_MaterialStateManagerInitializer** s,
    mgis_bv_MaterialDataManagerInitializer* const d) {
  if (d == nullptr) {
    return mgis_report_failure("invalid argument (material data manager is null)");
  }
  *s = &(d->s1);
  return mgis_report_success();
}  // end of mgis_bv_material_data_manager_get_state_1_initializer

mgis_status mgis_bv_free_material_data_manager_initializer(
    mgis_bv_MaterialDataManagerInitializer** d) {
  try {
    delete *d;
    *d = nullptr;
  } catch (...) {
    *d = nullptr;
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
} // end of mgis_bv_free_material_data_manager_initializer

mgis_status mgis_bv_create_material_data_manager(
    mgis_bv_MaterialDataManager** d,
    const mgis_bv_Behaviour* const b,
    const mgis_size_type n) {
  *d = nullptr;
  if (b == nullptr) {
    return mgis_report_failure(
        "mgis_bv_create_material_data_manager: "
        "null behaviour");
  }
  try {
    *d = new mgis::behaviour::MaterialDataManager(*b, n);
    if (*d == nullptr) {
      return mgis_report_failure(
          "mgis_bv_create_material_data_manager: "
          "memory allocation failed");
    }
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_create_material_data_manager

mgis_status mgis_bv_create_material_data_manager_with_initializer(
    mgis_bv_MaterialDataManager** d,
    const mgis_bv_Behaviour* const b,
    const mgis_size_type n,
    const mgis_bv_MaterialDataManagerInitializer* const i) {
  *d = nullptr;
  if (b == nullptr) {
    return mgis_report_failure(
        "mgis_bv_create_material_data_manager_with_initializer: "
        "null behaviour");
  }
  try {
    *d = new mgis::behaviour::MaterialDataManager(*b, n, *i);
    if (*d == nullptr) {
      return mgis_report_failure(
          "mgis_bv_create_material_data_manager_with_initializer: "
          "memory allocation failed");
    }
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_create_material_data_manager_with_initializer

mgis_status mgis_bv_material_data_manager_get_state_0(
    mgis_bv_MaterialStateManager** s, mgis_bv_MaterialDataManager* const d) {
  if (d == nullptr) {
    return mgis_report_failure("invalid argument (material data manager is null)");
  }
  *s = &(d->s0);
  return mgis_report_success();
}  // end of mgis_bv_material_data_manager_get_state_0

mgis_status mgis_bv_material_data_manager_get_state_1(
    mgis_bv_MaterialStateManager** s, mgis_bv_MaterialDataManager* const d) {
  if (d == nullptr) {
    return mgis_report_failure("invalid argument (material data manager is null)");
  }
  *s = &(d->s1);
  return mgis_report_success();
}  // end of mgis_bv_material_data_manager_get_state_1

mgis_status mgis_bv_material_data_manager_get_tangent_operator(
    mgis_real** const K, mgis_bv_MaterialDataManager* const d) {
  if (d == nullptr) {
    *K = nullptr;
    return mgis_report_failure("invalid argument (behaviour data is null)");
  }
  auto* const Kv = d->K.data();
  if (Kv == nullptr) {
    *K = nullptr;
    return mgis_report_failure("no tangent operator defined");
  }
  *K = &(Kv[0]);
  return mgis_report_success();
} // end of mgis_bv_material_data_manager_get_tangent_operator
  
mgis_status mgis_bv_update_material_data_manager(mgis_bv_MaterialDataManager* const d){
  try {
    mgis::behaviour::update(*d);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
} // end of mgis_bv_update_material_data_manager

mgis_status mgis_bv_revert_material_data_manager(mgis_bv_MaterialDataManager* const d){
  try {
    mgis::behaviour::revert(*d);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
} // end of mgis_bv_revert_material_data_manager

mgis_status mgis_bv_free_material_data_manager(mgis_bv_MaterialDataManager** d){
  try {
    delete *d;
    *d = nullptr;
  } catch (...) {
    *d = nullptr;
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
} // end of mgis_bv_free_material_data_manager

}  // end of extern "C"
