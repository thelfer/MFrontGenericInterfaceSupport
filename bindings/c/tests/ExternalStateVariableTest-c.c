/*!
 * \file   bindings/c/tests/ExternalStateVariableTest-c.c
 * \brief    
 * \author Thomas Helfer
 * \date   30/11/2021
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
#include "MGIS/Behaviour/Behaviour.h"
#include "MGIS/Behaviour/BehaviourData.h"
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

static void check_external_state_variable(
    const mgis_bv_Behaviour* const b,
    const char* const expected_variable_name,
    const mgis_bv_VariableType expected_variable_type,
    const mgis_size_type expected_variable_offset,
    const mgis_size_type i) {
  const char* name;
  mgis_bv_VariableType type;
  mgis_size_type offset;
  check_status(mgis_bv_behaviour_get_external_state_variable_name(&name, b, i));
  check(strcmp(name, expected_variable_name), "invalid external state variable name");
  check_status(mgis_bv_behaviour_get_external_state_variable_type(&type, b, i));
  check(type == expected_variable_type, "invalid external state variable type");
  check_status(
      mgis_bv_behaviour_get_external_state_variable_offset(&offset, b, i));
  check(offset == expected_variable_offset, "invalid external state variable offset");
}  // end of check_external_state_variable

static void check_behaviour(const mgis_bv_Behaviour* const b){
  // behaviour name
  const char* bn;
  check_status(mgis_bv_behaviour_get_behaviour_name(&bn, b));
  check_string(bn, "TensorialExternalStateVariableTest",
               "invalid behaviour name");
  // hypothesis
  const char* h;
  check_status(mgis_bv_behaviour_get_hypothesis(&h, b));
  check_string(h, "Tridimensional", "invalid hypothesis");
  // source
  const char* source;
  check_status(mgis_bv_behaviour_get_source(&source, b));
  check_string(source, "TensorialExternalStateVariableTest.mfront",
               "invalid source");
  // version
  const char* v;
  check_status(mgis_bv_behaviour_get_tfel_version(&v, b));
  check_string(v, TFEL_VERSION, "invalid TFEL version");
  // number of state variables
  mgis_size_type s;
  check_status(mgis_bv_behaviour_get_number_of_external_state_variables(&s, b));
  check(s == 10, "invalid number of external state variables");
  // type of the external state variables
  check_external_state_variable(b, "Temperature", MGIS_BV_SCALAR, 0, 0);
  check_external_state_variable(b, "v_esv", MGIS_BV_VECTOR, 1, 1);
  check_external_state_variable(b, "v2_esv[0]", MGIS_BV_VECTOR, 4, 2);
  check_external_state_variable(b, "v2_esv[1]", MGIS_BV_VECTOR, 7, 3);
}

static void check_behaviour_data(const mgis_bv_Behaviour* const b){
  mgis_bv_BehaviourData* d;
  mgis_bv_State *s1;
  mgis_bv_BehaviourDataView v;
  mgis_size_type i;
  mgis_real* isvs;
  mgis_real* esvs;
  int r;  // behaviour integration result
  const mgis_real eps = 1.e-14;
  //
  check_status(mgis_bv_allocate_behaviour_data(&d, b));
  // state at the end of the time step
  check_status(mgis_bv_behaviour_data_get_state_1(&s1, d));
  //
  const mgis_real v_esv_values[3] = {1, 2, 3};
  mgis_bv_state_set_external_state_variable_by_name(s1, "v_esv", v_esv_values);
  // initializing the view
  check_status(mgis_bv_make_behaviour_data_view(&v, d));
  check_status(mgis_bv_integrate(&r, &v, b));
  //
  mgis_bv_state_get_internal_state_variables(&isvs, s1);
  mgis_bv_state_get_external_state_variables(&esvs, s1);
  //
  for (i = 0; i != 3; ++i) {
    check(fabs(esvs[i + 1] - v_esv_values[i]) < eps,
          "invalid external state variable value");
    check(fabs(isvs[i] - v_esv_values[i]) < eps,
          "invalid internal state variable value");
  }
  // clean-up
  check_status(mgis_bv_free_behaviour_data(&d));
}

static void initialize_state(mgis_bv_MaterialStateManager* s,
                             mgis_real* const zeros) {
  check_status(
      mgis_bv_material_state_manager_set_uniform_scalar_material_property(
          s, "YoungModulus", 150e9));
  check_status(
      mgis_bv_material_state_manager_set_uniform_scalar_material_property(
          s, "PoissonRatio", 0.3));
  mgis_bv_material_state_manager_set_uniform_scalar_external_state_variable(
      s, "Temperature", 293.15);
  mgis_bv_material_state_manager_set_non_uniform_external_state_variable(
      s, "s_esv", zeros, MGIS_BV_EXTERNAL_STORAGE);
  mgis_bv_material_state_manager_set_non_uniform_external_state_variable(
      s, "s2_esv[0]", zeros, MGIS_BV_EXTERNAL_STORAGE);
  mgis_bv_material_state_manager_set_non_uniform_external_state_variable(
      s, "s2_esv[1]", zeros, MGIS_BV_EXTERNAL_STORAGE);
  mgis_bv_material_state_manager_set_non_uniform_external_state_variable(
      s, "t_esv", zeros, MGIS_BV_EXTERNAL_STORAGE);
  mgis_bv_material_state_manager_set_non_uniform_external_state_variable(
      s, "t2_esv[0]", zeros, MGIS_BV_EXTERNAL_STORAGE);
  mgis_bv_material_state_manager_set_non_uniform_external_state_variable(
      s, "t2_esv[1]", zeros, MGIS_BV_EXTERNAL_STORAGE);
}

static void check_material_data_manager(const mgis_bv_Behaviour* const b) {
  const mgis_real eps = 1.e-14;
  mgis_bv_MaterialDataManager* d;
  mgis_bv_MaterialStateManager *s0, *s1;
  mgis_real* isvs;
  mgis_real dt = 2;
  mgis_size_type i;
  int r;
  //
  check_status(mgis_bv_create_material_data_manager(&d, b, 2));
  check_status(mgis_bv_material_data_manager_get_state_0(&s0, d));
  check_status(mgis_bv_material_data_manager_get_state_1(&s1, d));
  //
  mgis_real zeros[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  mgis_real v_esv_values[3] = {1, 2, 3};
  mgis_real v2_esv0_values[9] = {1, 2, 3, 7, 6, 5};
  mgis_real v2_esv1_values[9] = {9, 8, 2, 3, 6, 4};
  //
  initialize_state(s0, zeros);
  initialize_state(s1, zeros);
  //
  mgis_bv_material_state_manager_get_internal_state_variables(&isvs, s1);
  // v_esv is uniform
  mgis_bv_material_state_manager_set_uniform_external_state_variable(
      s0, "v_esv", zeros, MGIS_BV_EXTERNAL_STORAGE);
  mgis_bv_material_state_manager_set_uniform_external_state_variable(
      s1, "v_esv", v_esv_values, MGIS_BV_LOCAL_STORAGE);
  // v2_esv[0] and v2_esv[1] are not uniform
  mgis_bv_material_state_manager_set_non_uniform_external_state_variable(
      s0, "v2_esv[0]", zeros, MGIS_BV_EXTERNAL_STORAGE);
  mgis_bv_material_state_manager_set_non_uniform_external_state_variable(
      s1, "v2_esv[0]", v2_esv0_values, MGIS_BV_EXTERNAL_STORAGE);
  mgis_bv_material_state_manager_set_non_uniform_external_state_variable(
      s0, "v2_esv[1]", zeros, MGIS_BV_EXTERNAL_STORAGE);
  mgis_bv_material_state_manager_set_non_uniform_external_state_variable(
      s1, "v2_esv[1]", v2_esv1_values, MGIS_BV_EXTERNAL_STORAGE);
  // checking if the external storage do work as expected
  v2_esv1_values[3] = -1;
  //
  check_status(mgis_bv_integrate_material_data_manager_part(
      &r, d, MGIS_BV_INTEGRATION_NO_TANGENT_OPERATOR, dt, 0, 2));
  //
  for (i = 0; i != 3; ++i) {
    check(fabs(isvs[i] - v_esv_values[i]) < eps,
          "invalid internal state variable value");
    check(fabs(isvs[54 + i] - v_esv_values[i]) < eps,
          "invalid internal state variable value");
    check(fabs(isvs[i + 3] - v2_esv0_values[i]) < eps,
          "invalid internal state variable value");
    check(fabs(isvs[57 + i] - v2_esv0_values[i + 3]) < eps,
          "invalid internal state variable value");
    check(fabs(isvs[i + 6] - v2_esv1_values[i]) < eps,
          "invalid internal state variable value");
    check(fabs(isvs[60 + i] - v2_esv1_values[i + 3]) < eps,
          "invalid internal state variable value");
  }
  // clean-up
  check_status(mgis_bv_free_material_data_manager(&d));
}  // end of check_material_data_manager

int main(const int argc, const char* const* argv) {
  mgis_bv_Behaviour* b;
  if (check(argc == 2, "expected two arguments") == 0) {
    return EXIT_FAILURE;
  }
  check_status(mgis_bv_load_behaviour(
      &b, argv[1], "TensorialExternalStateVariableTest", "Tridimensional"));
  //
  check_behaviour_data(b);
  check_material_data_manager(b);
  // clean-up
  check_status(mgis_bv_free_behaviour(&b));
  return test_status;
} // end of main
