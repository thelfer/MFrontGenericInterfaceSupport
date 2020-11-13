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
#include "MGIS/Behaviour/Behaviour.h"

int test_status = EXIT_SUCCESS;

static void check_status(const mgis_status s) {
  if (s.exit_status != MGIS_SUCCESS) {
    fprintf(stderr, "invalid function call: %s\n", s.msg);
    exit(EXIT_FAILURE);
  }
}  // end of check_status

static int check(const int b, const char* const e) {
  if (b == 0) {
    test_status = EXIT_FAILURE;
    fprintf(stderr, "%s\n", e);
  }
  return b;
}  // end of check

static int check_string(const char* const s1,
                        const char* const s2,
                        const char* const e) {
  const int b = strcmp(s1, s2) == 0;
  if (b == 0) {
    test_status = EXIT_FAILURE;
    fprintf(stderr, "%s (expected '%s', got '%s')\n", e, s2, s1);
  }
  return b;
}  // end of check_string

int main(const int argc, const char* const* argv) {
  if (!check(argc == 3, "expected three arguments")) {
    return EXIT_FAILURE;
  }
  int btype;
  check_status(
      mgis_bv_is_standard_finite_strain_behaviour(&btype, argv[1], argv[2]));
  check(btype, "invalid type");
  mgis_bv_FiniteStrainBehaviourOptions* o;
  check_status(mgis_bv_create_finite_strain_behaviour_options(&o));
  mgis_bv_Behaviour* b;
  check_status(mgis_bv_load_finite_strain_behaviour(&b, o, argv[1], argv[2],
                                                    "Tridimensional"));
  mgis_size_type mps_size;
  check_status(
      mgis_bv_behaviour_get_number_of_material_properties(&mps_size, b));
  if (check(mps_size == 16, "invalid number of material properties")) {
    const char* mps[16] = {"YoungModulus1",
                           "YoungModulus2",
                           "YoungModulus3",
                           "PoissonRatio12",
                           "PoissonRatio23",
                           "PoissonRatio13",
                           "ShearModulus12",
                           "ShearModulus23",
                           "ShearModulus13",
                           "m",
                           "K",
                           "C",
                           "R0",
                           "Q",
                           "b",
                           "d1"};
    for (mgis_size_type i = 0; i != 16; ++i) {
      const char* mp;
      check_status(mgis_bv_behaviour_get_material_property_name(&mp, b, i));
      check_string(mp, mps[i], "invalid name for the material property");
    }
  }
  mgis_bv_free_behaviour(&b);
  return test_status;
}  // end of main
