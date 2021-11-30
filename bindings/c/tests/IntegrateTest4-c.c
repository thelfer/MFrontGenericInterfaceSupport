/*!
 * \file   IntegrateTest4-c.c
 * \brief    
 * \author Thomas Helfer
 * \date   21/09/2018
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
#include <stdlib.h>
#include "MGIS/ThreadPool.h"
#include "MGIS/Model/Model.h"
#include "MGIS/Behaviour/MaterialDataManager.h"
#include "MGIS/Behaviour/Integrate.h"

int test_status = EXIT_SUCCESS;
mgis_model_Model* model = NULL;
mgis_ThreadPool* p = NULL;
mgis_bv_MaterialDataManager* m = NULL;

static void check_status(const mgis_status s) {
  if (s.exit_status != MGIS_SUCCESS) {
    fprintf(stderr, "invalid function call: %s\n", s.msg);
    mgis_model_free_model(&model);
    mgis_bv_free_material_data_manager(&m);
    mgis_free_thread_pool(&p);
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

int main(const int argc, const char* const* argv) {
  if (check(argc == 2, "expected two arguments") == 0) {
    return EXIT_FAILURE;
  }
  const mgis_real eps = 1e-10;
  const mgis_real dt =0.1;
  mgis_size_type o;
  mgis_real* isvs0; /* internal state variables at the beginning of the time step */
  mgis_real* isvs1; /* internal state variables at the end of the time step */
  mgis_size_type isvs_stride; /* internal state variables stride */
  const mgis_size_type n = 100;
  mgis_bv_MaterialStateManager*
      s0; /* state at the beginning of the time step */
  mgis_bv_MaterialStateManager* s1; /* state at the end of the time step */
  mgis_real xi[11];      /* values of 'x' for the first integration point */
  mgis_real xe[11];      /* values of 'x' for the last integration point */
  mgis_size_type ni,ne;
  mgis_size_type idx;
  mgis_size_type i;
  mgis_real A;
  mgis_real t;
  mgis_real x_ref;
  int r;
  check_status(mgis_create_thread_pool(&p, 2));
  check_status(mgis_model_load(&model, argv[1], "ode_rk54", "Tridimensional"));
  check_status(mgis_bv_behaviour_get_parameter_default_value(&A, model, "A"));
  check_status(mgis_bv_create_material_data_manager(&m, model, 100));
  check_status(mgis_bv_behaviour_get_internal_state_variable_offset_by_name(
      &o, model, "x"));
  check_status(mgis_bv_material_data_manager_get_state_0(&s0, m));
  check_status(mgis_bv_material_data_manager_get_state_1(&s1, m));
  /* initialize the internal state variable */
  check_status(
      mgis_bv_material_state_manager_get_internal_state_variables(&isvs0, s0));
  check_status(
      mgis_bv_material_state_manager_get_internal_state_variables(&isvs1, s1));
  check_status(
      mgis_bv_material_state_manager_get_internal_state_variables_stride(
          &isvs_stride, s1));
  for (idx = 0; idx != n; ++idx) {
    isvs1[idx * isvs_stride + o] = 1;
  }
  /* initialize the external state variable */
  check_status(
      mgis_bv_material_state_manager_set_uniform_scalar_external_state_variable(
          s1, "Temperature", 293.15));
  /* copy s1 in s0 */
  check_status(mgis_bv_update_material_data_manager(m));
  // integration */
  ni = o;
  ne = (n - 1) * isvs_stride + o;
  xi[0] = isvs0[ni];
  xe[0] = isvs0[ne];
  for (i = 0; i != 10; ++i) {
    check_status(mgis_bv_integrate_material_data_manager(
        &r, p, m, MGIS_BV_INTEGRATION_NO_TANGENT_OPERATOR, dt));
    check_status(mgis_bv_update_material_data_manager(m));
    xi[i + 1] = isvs1[ni];
    xe[i + 1] = isvs1[ne];
  }
  check_status(mgis_model_free_model(&model));
  check_status(mgis_bv_free_material_data_manager(&m));
  check_status(mgis_free_thread_pool(&p));
  t = 0;
  for (i = 0; i != 11; ++i) {
    x_ref = exp(-A * t);
    if (fabs(xi[i] - x_ref) > eps) {
      fprintf(stderr,
              "IntegrateTest: invalid value for x "
              "at the first integration point"
              "(expected '%g', computed '%g')\n",
              x_ref, xi[i]);
      return EXIT_FAILURE;
    } 
    if (fabs(xe[i] - x_ref) > eps) {
      fprintf(stderr,
              "IntegrateTest: invalid value for x "
              "at the first integration point"
              "(expected '%g', computed '%g')\n",
              x_ref, xe[i]);
      return EXIT_FAILURE;
    }
    t += dt;
  }
  return EXIT_SUCCESS;
}
