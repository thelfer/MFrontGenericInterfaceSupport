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
#include <string.h>
#include <stdlib.h>
#include "MGIS/Behaviour/Description.h"

int test_status = EXIT_SUCCESS;

static void check_status(const mgis_status s) {
  if (s.exit_status != MGIS_SUCCESS) {
    fprintf(stderr, "invalid function call: %s\n",s.msg);
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
  if (check(argc == 3, "expected three arguments") == 0) {
    return EXIT_FAILURE;
  }
  mgis_bv_Description *d;
  check_status(
      mgis_bv_load_description(&d, argv[1], argv[2], "Tridimensional"));
  // behaviour name
  const char* b;
  check_status(mgis_bv_get_behaviour_name(&b,d));
  check_string(b, "gurson", "invalid behaviour name");
  // hypothesis
  const char* h;
  check_status(mgis_bv_get_hypothesis(&h,d));
  check_string(h, "Tridimensional", "invalid hypothesis");
  // source
  const char* s;
  check_status(mgis_bv_get_source(&s,d));
  check_string(s, "Gurson.mfront", "invalid source");
  // version
  const char* v;
  check_status(mgis_bv_get_tfel_version(&v,d));
  check_string(v, "3.2.0-dev", "invalid TFEL version");
  // material properties
  mgis_size_type mps_size;
  check_status(mgis_bv_get_number_of_material_properties(&mps_size, d));
  check(mps_size == 0, "invalid number of material properties");
  // internal state variables
  mgis_size_type isvs_size;
  check_status(mgis_bv_get_number_of_internal_state_variables(&isvs_size, d));
  if (check(isvs_size == 4, "invalid number of internal state variables")) {
    const char* eel;
    check_status(mgis_bv_get_internal_state_variable_name(&eel, d, 0));
    check_string(eel, "ElasticStrain", "invalid internal state variable name");
    mgis_bv_VariableType eel_t;
    check_status(mgis_bv_get_internal_state_variable_type(&eel_t, d, 0));
    check(eel_t == MGIS_BV_STENSOR,
          "invalid type for internal state variable 'ElasticStrain'");
    const char* p;
    check_status(mgis_bv_get_internal_state_variable_name(&p, d, 1));
    check_string(p, "EquivalentPlasticStrain",
                 "invalid internal state variable name");
    mgis_bv_VariableType p_t;
    check_status(mgis_bv_get_internal_state_variable_type(&p_t, d, 1));
    check(p_t == MGIS_BV_SCALAR,
          "invalid type for internal state variable 'EquivalentPlasticStrain'");
    const char* pm;
    check_status(mgis_bv_get_internal_state_variable_name(&pm, d, 2));
    check_string(pm, "MatrixEquivalentPlasticStrain",
                 "invalid internal state variable name");
    mgis_bv_VariableType pm_t;
    check_status(mgis_bv_get_internal_state_variable_type(&pm_t, d, 2));
    check(pm_t == MGIS_BV_SCALAR,
          "invalid type for internal state variable 'MatrixEquivalentPlasticStrain'");
    const char* f;
    check_status(mgis_bv_get_internal_state_variable_name(&f, d, 3));
    check_string(f, "Porosity", "invalid internal state variable name");
    mgis_bv_VariableType f_t;
    check_status(mgis_bv_get_internal_state_variable_type(&f_t, d, 2));
    check(f_t == MGIS_BV_SCALAR,
          "invalid type for internal state variable 'Porosity'");
  }
  // external state variables
  mgis_size_type esvs_size;
  check_status(mgis_bv_get_number_of_external_state_variables(&esvs_size, d));
  if (check(esvs_size == 1, "invalid number of external state variables")) {
    const char* T;
    check_status(mgis_bv_get_external_state_variable_name(&T, d, 0));
    check_string(T, "Temperature", "invalid external state variable name");
  }
  mgis_bv_free_description(&d);
  return test_status;
}  // end of main
