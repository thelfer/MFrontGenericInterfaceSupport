/*!
 * \file   bindings/c/src/BehaviourDataView.cxx
 * \brief    
 * \author th202608
 * \date   02/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include "MGIS/Behaviour/BehaviourDataView.h"

extern "C" {

mgis_status mgis_bv_make_behaviour_data_view(mgis_bv_BehaviourDataView* const v,
                                             mgis_bv_BehaviourData* const d) {
  if(v==nullptr){
    return mgis_report_failure("uninitialized view");
  }
  if(d==nullptr){
    return mgis_report_failure("uninitialized behaviour data");
  }
  try {
    *v = make_view(*d);
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_make_behaviour_data_view

} // end of extern "C"
