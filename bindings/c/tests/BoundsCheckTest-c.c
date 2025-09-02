/*!
 * \file   BoundsCheckTest-c.c
 * \brief    
 * \author Thomas Helfer
 * \date   13/11/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "MGIS/Behaviour/Integrate.h"

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
  typedef long double lreal;
  const lreal yg_min = 100.e9;
  const lreal yg_max = 200.e9;
  const lreal eps = 1.e-14;
  mgis_bv_Behaviour* b;
  if (check(argc == 2, "expected two arguments") == 0) {
    return EXIT_FAILURE;
  }
  check_status(
      mgis_bv_load_behaviour(&b, argv[1], "BoundsCheckTest", "Tridimensional"));
  // behaviour name
  const char* bn;
  check_status(mgis_bv_behaviour_get_behaviour_name(&bn, b));
  check_string(bn, "BoundsCheckTest", "invalid behaviour name");
  // hypothesis
  const char* h;
  check_status(mgis_bv_behaviour_get_hypothesis(&h, b));
  check_string(h, "Tridimensional", "invalid hypothesis");
  // source
  const char* s;
  check_status(mgis_bv_behaviour_get_source(&s, b));
  check_string(s, "BoundsCheckTest.mfront", "invalid source");
  // version
  const char* v;
  check_status(mgis_bv_behaviour_get_tfel_version(&v, b));
  check_string(v, TFEL_VERSION, "invalid TFEL version");
  // bounds
  int hasBounds;
  check_status(mgis_bv_behaviour_has_bounds(&hasBounds, b, "YoungModulus"));
  check(hasBounds, "'YoungModulus' shall have bounds");
  int hasLowerBound;
  check_status(mgis_bv_behaviour_has_lower_bound(&hasLowerBound, b, "YoungModulus"));
  check(hasLowerBound, "'YoungModulus' shall have a lower bound");
  long double lowerBound;
  check_status(
      mgis_bv_behaviour_get_lower_bound(&lowerBound, b, "YoungModulus"));
  check(fabsl(lowerBound - yg_min) < eps * yg_min,
        "invalid lower bound for 'YoungModulus'");
  int hasUpperBound;
  check_status(mgis_bv_behaviour_has_upper_bound(&hasUpperBound, b, "YoungModulus"));
  check(hasUpperBound, "'YoungModulus' shall have an upper bound");
  long double upperBound;
  check_status(
      mgis_bv_behaviour_get_upper_bound(&upperBound, b, "YoungModulus"));
  check(fabsl(upperBound - yg_max) < eps * yg_max,
        "invalid upper bound for 'YoungModulus'");
  // physical bounds
  int hasPhysicalBounds;
  check_status(mgis_bv_behaviour_has_physical_bounds(&hasPhysicalBounds, b, "YoungModulus"));
  check(hasPhysicalBounds, "'YoungModulus' shall have bounds");
  int hasLowerPhysicalBound;
  check_status(mgis_bv_behaviour_has_lower_physical_bound(&hasLowerPhysicalBound, b, "YoungModulus"));
  check(hasLowerPhysicalBound, "'YoungModulus' shall have a lower bound");
  long double lowerPhysicalBound;
  check_status(mgis_bv_behaviour_get_lower_physical_bound(&lowerPhysicalBound,
                                                          b, "YoungModulus"));
  check(fabsl(lowerPhysicalBound) < eps * yg_min,
        "invalid physical lower bound for 'YoungModulus'");
  int hasUpperPhysicalBound;
  check_status(mgis_bv_behaviour_has_upper_physical_bound(&hasUpperPhysicalBound, b, "YoungModulus"));
  check(!hasUpperPhysicalBound, "'YoungModulus' shall have an upper bound");
  // clean-up
  check_status(mgis_bv_free_behaviour(&b));
  return test_status;
} // end of main