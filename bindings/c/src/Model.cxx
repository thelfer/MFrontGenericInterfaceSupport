/*!
 * \file   Model.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/10/2021
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include "MGIS/Model/Model.h"

extern "C" {

mgis_status mgis_model_load(mgis_model_Model** ptr,
                            const char* const l,
                            const char* const m,
                            const char* const h) {
  *ptr = nullptr;
  try {
    const auto model = mgis::model::load(l, m, mgis::behaviour::fromString(h));
    *ptr = new mgis::model::Model(std::move(model));
    if (*ptr == nullptr) {
      return mgis_report_failure(
          "mgis_model_load: "
          "memory allocation failed");
    }
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_model_load

mgis_status mgis_model_free_model(mgis_model_Model** m) {
  return mgis_bv_free_behaviour(m);
}  // end of mgis_model_free_model

}  // end of extern "C"
