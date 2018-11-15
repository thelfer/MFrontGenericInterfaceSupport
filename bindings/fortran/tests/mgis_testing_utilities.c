/*!
 * \file   mgis_testing_utilities.c
 * \brief    
 * \date   15/11/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include<stdlib.h>

const char* mgis_get_mfront_behaviour_test_library_path(){
  return getenv("MGIS_TEST_BEHAVIOURS_LIBRARY");
} // end of mgis_get_mfront_behaviour_test_library_path

const char* mgis_get_tfel_version(){
  return getenv("MGIS_TEST_TFEL_VERSION");
} // end of mgis_get_mfront_behaviour_test_library
