/*!
 * \file   MFrontGenericBehaviourInterfaceTest-c.c
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

#include <cstdlib>
#include "MGIS/Behaviour/Description.h"

int test_status = EXIT_SUCCESS;

int check(const int b, const char* const e) {
  if (!b) {
    fprintf(stderr, "%s\n", e);
  }
  return b;
}  // end of check

int main(const int argc, const char* const* argv) {
  if (check(argc != 3, "expected three arguments")) {
    return EXIT_FAILURE;
  }
  mgis_behaviour_Description *d;
  mgis_behaviour_load(&d, argv[1], argv[2], "Tridimensional");
 return test_status;
}  // end of main