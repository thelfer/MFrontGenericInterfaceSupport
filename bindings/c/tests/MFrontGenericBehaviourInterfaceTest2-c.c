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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MGIS/Behaviour/Description.h"

int test_status = EXIT_SUCCESS;

static void check_status(const mgis_status s) {
  if (s.exit_status != MGIS_SUCCESS) {
    fprintf(stderr, "invalid function call\n");
    exit(EXIT_FAILURE);
  }
} // end of check_status

static int check(const int b, const char* const e) {
  if (b == 0) {
    test_status = EXIT_FAILURE;
    fprintf(stderr, "%s\n", e);
  }
  return b;
}  // end of check

int main(const int argc, const char* const* argv) {
  if (!check(argc == 3, "expected three arguments")) {
    return EXIT_FAILURE;
  }
  mgis_behaviour_Description *d;
  check_status(
      mgis_behaviour_load_description(&d, argv[1], argv[2], "Tridimensional"));
  mgis_size_type mps_size;
  check_status(mgis_behaviour_get_number_of_material_properties(&mps_size, d));
  if (check(mps_size == 16, "invalid number of material properties")) {
    const char *mp1;
    check_status(mgis_behaviour_get_material_property_name(&mp1, d, 0));
    check(strcmp(mp1, "YoungModulus1") == 0,
          "invalid name for the material property 'YoungModulus1'");
  }
  mgis_behaviour_free_description(&d);
  return test_status;
}  // end of main
