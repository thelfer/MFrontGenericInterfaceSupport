/*!
 * \file   InitializeFunctionAndPostProcessingTest-c.c
 * \brief  tests of the initialize functions and of the post-processings
 * \date   27/10/2026
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "MGIS/Behaviour/State.h"
#include "MGIS/Behaviour/BehaviourData.h"
#include "MGIS/Behaviour/MaterialDataManager.h"
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

static void test_initialize_function(const char* const l) {
  const char* const f = "StressFromInitialPressure";
  const mgis_real pr = -1.2e5;
  const mgis_real eps = -10 * pr * DBL_EPSILON;
  mgis_bv_Behaviour* b;
  mgis_bv_BehaviourData* d;
  mgis_bv_MaterialDataManager* m;
  mgis_bv_State *s0, *s1;
  mgis_bv_MaterialStateManager *ms0, *ms1;
  mgis_bv_VariableType t;
  const char* n;
  mgis_size_type c, i, j, stride;
  mgis_real* sig;
  int r;
  check_status(mgis_bv_load_behaviour(&b, l, "InitializeFunctionTest",
                                      "Tridimensional"));
  // description
  check_status(mgis_bv_behaviour_get_number_of_initialize_functions(&c, b));
  check(c == 3, "invalid number of initialize functions");
  check_status(mgis_bv_behaviour_get_initialize_function_name(&n, b, 2));
  check(strcmp(n, f) == 0, "invalid initialize function name");
  check_status(mgis_bv_behaviour_get_initialize_function_inputs_size(&c, b, f));
  check(c == 1, "invalid size of the inputs");
  check_status(
      mgis_bv_behaviour_get_number_of_initialize_function_inputs(&c, b, f));
  check(c == 1, "invalid number of inputs");
  check_status(
      mgis_bv_behaviour_get_initialize_function_input_name(&n, b, f, 0));
  check(strcmp(n, "pr") == 0, "invalid input name");
  check_status(
      mgis_bv_behaviour_get_initialize_function_input_type(&t, b, f, 0));
  check(t == MGIS_BV_SCALAR, "invalid input type");
  check_status(
      mgis_bv_behaviour_get_initialize_function_input_offset(&c, b, f, 0));
  check(c == 0, "invalid input offset");
  // execution on a behaviour data
  check_status(mgis_bv_allocate_behaviour_data(&d, b));
  check_status(mgis_bv_behaviour_data_get_state_0(&s0, d));
  check_status(mgis_bv_behaviour_data_get_state_1(&s1, d));
  check_status(mgis_bv_state_set_scalar_external_state_variable_by_name(
      s0, "Temperature", 293.15));
  check_status(mgis_bv_state_set_scalar_external_state_variable_by_name(
      s1, "Temperature", 293.15));
  check_status(mgis_bv_execute_initialize_function_2(&r, d, b, f, &pr, 1));
  check(r == 0, "invalid result of the initialize function");
  check_status(
      mgis_bv_state_get_thermodynamic_force_by_name(&sig, s1, "Stress"));
  for (i = 0; i != 6; ++i) {
    check(fabs(sig[i] - (i < 3 ? pr : 0)) < eps, "invalid stress value");
  }
  check(mgis_bv_execute_initialize_function_2(&r, d, b, "Failure", NULL, 0)
                .exit_status == MGIS_FAILURE,
        "the failure of the initialize function is not reported");
  // execution on a material data manager, with uniform inputs
  check_status(mgis_bv_create_material_data_manager(&m, b, 2));
  check_status(mgis_bv_material_data_manager_get_state_0(&ms0, m));
  check_status(mgis_bv_material_data_manager_get_state_1(&ms1, m));
  check_status(
      mgis_bv_material_state_manager_set_uniform_scalar_external_state_variable(
          ms0, "Temperature", 293.15));
  check_status(
      mgis_bv_material_state_manager_set_uniform_scalar_external_state_variable(
          ms1, "Temperature", 293.15));
  check_status(mgis_bv_execute_initialize_function_material_data_manager_part(
      &r, m, f, &pr, 1, 0, 2));
  check(r == 1, "invalid result of the initialize function");
  check_status(
      mgis_bv_material_state_manager_get_thermodynamic_forces(&sig, ms1));
  check_status(mgis_bv_material_state_manager_get_thermodynamic_forces_stride(
      &stride, ms1));
  for (j = 0; j != 2; ++j) {
    for (i = 0; i != 6; ++i) {
      check(fabs(sig[j * stride + i] - (i < 3 ? pr : 0)) < eps,
            "invalid stress value");
    }
  }
  check_status(mgis_bv_free_material_data_manager(&m));
  check_status(mgis_bv_free_behaviour_data(&d));
  check_status(mgis_bv_free_behaviour(&b));
}  // end of test_initialize_function

static void test_post_processing(const char* const l) {
  const char* const p = "PrincipalStrain";
  const mgis_real e[6] = {1.3e-2, 1.2e-2, 1.4e-2, 0., 0., 0.};
  const mgis_real e2[3] = {1.2e-2, 1.3e-2, 1.4e-2};
  const mgis_real eps = 10 * DBL_EPSILON;
  mgis_bv_Behaviour* b;
  mgis_bv_BehaviourData* d;
  mgis_bv_BehaviourDataView v;
  mgis_bv_MaterialDataManager* m;
  mgis_ThreadPool* tp;
  mgis_bv_State *s0, *s1;
  mgis_bv_MaterialStateManager* ms1;
  mgis_bv_VariableType t;
  const char* n;
  mgis_size_type c, i, j, stride;
  mgis_real* g;
  mgis_real outputs[6];
  int r;
  check_status(
      mgis_bv_load_behaviour(&b, l, "PostProcessingTest", "Tridimensional"));
  // description
  check_status(mgis_bv_behaviour_get_number_of_post_processings(&c, b));
  check(c == 1, "invalid number of post-processings");
  check_status(mgis_bv_behaviour_get_post_processing_name(&n, b, 0));
  check(strcmp(n, p) == 0, "invalid post-processing name");
  check_status(mgis_bv_behaviour_get_post_processing_outputs_size(&c, b, p));
  check(c == 3, "invalid size of the outputs");
  check_status(
      mgis_bv_behaviour_get_number_of_post_processing_outputs(&c, b, p));
  check(c == 1, "invalid number of outputs");
  check_status(mgis_bv_behaviour_get_post_processing_output_name(&n, b, p, 0));
  check(strcmp(n, "PrincipalStrain") == 0, "invalid output name");
  check_status(mgis_bv_behaviour_get_post_processing_output_type(&t, b, p, 0));
  check(t == MGIS_BV_VECTOR_3D, "invalid output type");
  check_status(
      mgis_bv_behaviour_get_post_processing_output_offset(&c, b, p, 0));
  check(c == 0, "invalid output offset");
  // execution on a behaviour data view
  check_status(mgis_bv_allocate_behaviour_data(&d, b));
  check_status(mgis_bv_behaviour_data_get_state_0(&s0, d));
  check_status(mgis_bv_behaviour_data_get_state_1(&s1, d));
  check_status(mgis_bv_state_set_scalar_external_state_variable_by_name(
      s0, "Temperature", 293.15));
  check_status(mgis_bv_state_set_scalar_external_state_variable_by_name(
      s1, "Temperature", 293.15));
  check_status(mgis_bv_state_set_gradient_by_name(s1, "Strain", e));
  check_status(mgis_bv_make_behaviour_data_view(&v, d));
  check_status(mgis_bv_execute_post_processing(&r, outputs, 3, &v, b, p));
  check(r == 0, "invalid result of the post-processing");
  for (i = 0; i != 3; ++i) {
    check(fabs(outputs[i] - e2[i]) < eps, "invalid output value");
  }
  // execution on a material data manager, using a thread pool
  check_status(mgis_bv_create_material_data_manager(&m, b, 2));
  check_status(mgis_bv_material_data_manager_get_state_1(&ms1, m));
  check_status(
      mgis_bv_material_state_manager_set_uniform_scalar_material_property(
          ms1, "YoungModulus", 150e9));
  check_status(
      mgis_bv_material_state_manager_set_uniform_scalar_material_property(
          ms1, "PoissonRatio", 0.3));
  check_status(
      mgis_bv_material_state_manager_set_uniform_scalar_external_state_variable(
          ms1, "Temperature", 293.15));
  check_status(mgis_bv_update_material_data_manager(m));
  check_status(mgis_bv_material_state_manager_get_gradients(&g, ms1));
  check_status(
      mgis_bv_material_state_manager_get_gradients_stride(&stride, ms1));
  for (j = 0; j != 2; ++j) {
    for (i = 0; i != 6; ++i) {
      g[j * stride + i] = e[i];
    }
  }
  check_status(mgis_create_thread_pool(&tp, 2));
  check_status(mgis_bv_execute_post_processing_material_data_manager(
      &r, outputs, 6, tp, m, p));
  check(r == 1, "invalid result of the post-processing");
  for (j = 0; j != 2; ++j) {
    for (i = 0; i != 3; ++i) {
      check(fabs(outputs[3 * j + i] - e2[i]) < eps, "invalid output value");
    }
  }
  check_status(mgis_free_thread_pool(&tp));
  check_status(mgis_bv_free_material_data_manager(&m));
  check_status(mgis_bv_free_behaviour_data(&d));
  check_status(mgis_bv_free_behaviour(&b));
}  // end of test_post_processing

int main(const int argc, const char* const* argv) {
  if (check(argc == 2, "expected two arguments") == 0) {
    return EXIT_FAILURE;
  }
  test_initialize_function(argv[1]);
  test_post_processing(argv[1]);
  return test_status;
}  // end of main
