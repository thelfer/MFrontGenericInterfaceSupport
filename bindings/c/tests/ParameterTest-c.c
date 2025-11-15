/*!
 * \file   ParameterTest-c.c
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
  const mgis_real yg = 150e9;
  const mgis_real nu = 0.3;
  const mgis_real eps = 1.e-14;
  mgis_bv_Behaviour* b;
  if (check(argc == 2, "expected two arguments") == 0) {
    return EXIT_FAILURE;
  }
  check_status(
      mgis_bv_load_behaviour(&b, argv[1], "ParameterTest", "Tridimensional"));
  // behaviour name
  const char* bn;
  check_status(mgis_bv_behaviour_get_behaviour_name(&bn, b));
  check_string(bn, "ParameterTest", "invalid behaviour name");
  // hypothesis
  const char* h;
  check_status(mgis_bv_behaviour_get_hypothesis(&h, b));
  check_string(h, "Tridimensional", "invalid hypothesis");
  // source
  const char* s;
  check_status(mgis_bv_behaviour_get_source(&s, b));
  check_string(s, "ParameterTest.mfront", "invalid source");
  // version
  const char* tfel_version;
  check_status(mgis_bv_behaviour_get_tfel_version(&tfel_version, b));
  check_string(tfel_version, TFEL_VERSION, "invalid TFEL version");
  // number of real parameter
  mgis_size_type nparams;
  check_status(mgis_bv_behaviour_get_number_of_parameters(&nparams, b));
  if (check(nparams == 6u, "invalid number of parameters")) {
    const char* n;
    mgis_real v;
    check_status(mgis_bv_behaviour_get_parameter_name(&n, b, 0));
    check(strcmp(n, "YoungModulus") == 0, "invalid first parameter");
    check_status(mgis_bv_behaviour_get_parameter_default_value(&v, b, n));
    check(fabs(v-yg)<eps*yg,"invalid 'YoungModulus' default value");
    check_status(mgis_bv_behaviour_get_parameter_name(&n, b, 1));
    check(strcmp(n, "PoissonRatio") == 0, "invalid second parameter");
    check_status(mgis_bv_behaviour_get_parameter_default_value(&v, b, n));
    check(fabs(v-nu)<eps*nu,"invalid 'PoissonRatio' default value");
    check_status(mgis_bv_behaviour_get_parameter_name(&n, b, 2));
    check(strcmp(n, "ParametersArray[0]") == 0, "invalid third parameter");
    check_status(mgis_bv_behaviour_get_parameter_name(&n, b, 3));
    check(strcmp(n, "ParametersArray[1]") == 0, "invalid fourth parameter");
    check_status(mgis_bv_behaviour_get_parameter_name(&n, b, 4));
    check(strcmp(n, "minimal_time_step_scaling_factor") == 0, "invalid fifth parameter");
    check_status(mgis_bv_behaviour_get_parameter_name(&n, b, 5));
    check(strcmp(n, "maximal_time_step_scaling_factor") == 0, "invalid sixth parameter");
  }
  // clean-up
  check_status(mgis_bv_free_behaviour(&b));
  return test_status;
} // end of main
