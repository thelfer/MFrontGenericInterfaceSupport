/*!
 * \file   Description.cxx
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

#include <utility>
#include <stdexcept>
#include "MGIS/Behaviour/Description.h"

extern "C" {

mgis_status load(mgis_behaviour_Description** ptr,
                 const char* const l,
                 const char* const b,
                 const char* const h) {
  *ptr = nullptr;
  try {
    const auto d = mgis::behaviour::load(l, b, mgis::behaviour::fromString(h));
    *ptr = new mgis::behaviour::Description(std::move(d));
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
} // end of load

void mgis_behaviour_free_Description(mgis_behaviour_Description* const d) {
  std::free(d);
} // end of mgis_behaviour_free_Description

} // end of extern "C"