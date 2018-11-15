/*!
 * \file   Status.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   24/06/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <cstring>
#include <stdexcept>
#include <algorithm>
#include "MGIS/Status.h"

extern "C" {

mgis_status mgis_report_success() {
  return { MGIS_SUCCESS, nullptr };
}  // end of mgis_status mgis_report_success

mgis_status mgis_report_failure(const char* const e) {
  static thread_local char msg[512];
  ::strncpy(msg, e, 511);
  msg[511] = '\0';
  return {MGIS_FAILURE, msg};
} // end of mgis_status mgis_report_failure

mgis_status mgis_handle_cxx_exception() {
  try {
    throw;
  } catch (std::exception& e) {
    return mgis_report_failure(e.what());
  } catch (...) {
    return mgis_report_failure(
        "mgis_handle_cxx_exception: "
        "unknown exception");
  }
}  // end of mgis_handle_cxx_exception

}  // end of extern "C"
