/*!
 * \file   ThreadPool.cxx
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

#include "MGIS/ThreadPool.h"

mgis_status mgis_bv_create_thread_pool(mgis_ThreadPool** p,
                                       const mgis_size_type n) {
  *p = nullptr;
  try {
    *p = new mgis::ThreadPool(n);
    if (*p == nullptr) {
      return mgis_report_failure(
          "mgis_bv_create_thread_pool: "
          "memory allocation failed");
    }
  } catch (...) {
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
}  // end of mgis_bv_create_thread_pool

mgis_status mgis_bv_free_thread_pool(mgis_ThreadPool** p){
  try {
    delete *p;
    *p = nullptr;
  } catch (...) {
    *p = nullptr;
    return mgis_handle_cxx_exception();
  }
  return mgis_report_success();
} // end of mgis_bv_free_thread_pool

