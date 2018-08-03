/*!
 * \file   bindings/c/src/BehaviourData.cxx
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

#include "MGIS/Behaviour/BehaviourData.hxx"
#include "MGIS/Behaviour/BehaviourData.h"

extern "C" {

mgis_status mgis_bv_allocate_behaviour_data(mgis_bv_BehaviourData** d,
                                            const mgis_bv_Behaviour* const b) {
  if (b == nullptr) {
    *d = nullptr;
    return mgis_report_failure("invalid argument (behaviour is null)");
  }
  try {
    *d = new mgis::behaviour::BehaviourData(*b);
  } catch (...) {
    *d = nullptr;
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_allocate_behaviour_data

mgis_status mgis_bv_update_behaviour_data(mgis_bv_BehaviourData* const d){
  try {
    revert(*d);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
} // end of mgis_bv_update_behaviour_data

mgis_status mgis_bv_revert_behaviour_data(mgis_bv_BehaviourData* const d) {
  try {
    revert(*d);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
} // end of mgis_bv_revert_behaviour_data

mgis_status mgis_bv_behaviour_data_get_state_0(
    mgis_bv_State** s, mgis_bv_BehaviourData* const d) {
  *s = &(d->s0);
  return mgis_report_success();
}  // end of mgis_bv_behaviour_data_get_state_0

mgis_status mgis_bv_behaviour_data_get_state_1(
    mgis_bv_State** s, mgis_bv_BehaviourData* const d) {
  *s = &(d->s1);
  return mgis_report_success();
}  // end of mgis_bv_behaviour_data_get_state_1

mgis_status mgis_bv_free_behaviour_data(mgis_bv_BehaviourData** d) {
  try {
    delete *d;
   } catch (...) {
     return mgis_handle_cxx_exception();
  }
  *d = nullptr;
  return mgis_report_success();
}  // end of mgis_bv_free_behaviour_data

}  // end of extern "C"