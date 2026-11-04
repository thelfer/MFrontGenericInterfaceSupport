/*!
 * \file   bindings/c/tests/TangentOperatorTest-c.c
 * \brief  tests of the tangent and prediction operators returned by the
 * integrations of a material data manager and of a behaviour data, with the
 * `TangentKind` behaviour whose operator identifies the type of integration
 * \date   03/11/2026
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
#include "MGIS/Behaviour/Behaviour.h"
#include "MGIS/Behaviour/BehaviourData.h"
#include "MGIS/Behaviour/BehaviourDataView.h"
#include "MGIS/Behaviour/MaterialDataManager.h"
#include "MGIS/Behaviour/MaterialStateManager.h"
#include "MGIS/Behaviour/Integrate.h"

int test_status = EXIT_SUCCESS;
mgis_bv_Behaviour* b = NULL;
mgis_ThreadPool* p = NULL;
mgis_bv_MaterialDataManager* m = NULL;
mgis_bv_BehaviourData* d = NULL;

static void check_status(const mgis_status s) {
  if (s.exit_status != MGIS_SUCCESS) {
    fprintf(stderr, "invalid function call: %s\n", s.msg);
    mgis_bv_free_behaviour(&b);
    mgis_bv_free_material_data_manager(&m);
    mgis_bv_free_behaviour_data(&d);
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

/*
 * component (i, j) of the operator computed by the `TangentKind` behaviour
 * for the type of integration t and the value k of the material property
 */
static mgis_real expected_operator(const int t,
                                   const mgis_real k,
                                   const mgis_size_type i,
                                   const mgis_size_type j) {
  mgis_real v = (i == j) ? t * k : 0;
  if ((t > 0) && (i == 0) && (j == 1)) {
    v += k / 2;
  }
  return v;
}  // end of expected_operator

int main(const int argc, const char* const* argv) {
  const mgis_bv_IntegrationType types[7] = {
      MGIS_BV_PREDICTION_TANGENT_OPERATOR,
      MGIS_BV_PREDICTION_SECANT_OPERATOR,
      MGIS_BV_PREDICTION_ELASTIC_OPERATOR,
      MGIS_BV_INTEGRATION_ELASTIC_OPERATOR,
      MGIS_BV_INTEGRATION_SECANT_OPERATOR,
      MGIS_BV_INTEGRATION_TANGENT_OPERATOR,
      MGIS_BV_INTEGRATION_CONSISTENT_TANGENT_OPERATOR};
  mgis_real k[4] = {10, 20, 30, 40}; /* material property at each point */
  const mgis_size_type n = 4;        /* number of integration points */
  /* states at the beginning and at the end of the time step */
  mgis_bv_MaterialStateManager* s0;
  mgis_bv_MaterialStateManager* s1;
  /* states of the behaviour data */
  mgis_bv_State* ds0;
  mgis_bv_State* ds1;
  mgis_bv_BehaviourDataView v;
  mgis_real* g;    /* gradients */
  mgis_real* isvs; /* internal state variables */
  mgis_real* K;    /* tangent operators */
  mgis_size_type g_stride;
  mgis_size_type K_stride;
  mgis_size_type it, mode, idx, i, j;
  int r;
  if (check(argc == 2, "expected two arguments") == 0) {
    return EXIT_FAILURE;
  }
  check_status(mgis_create_thread_pool(&p, 2));
  check_status(
      mgis_bv_load_behaviour(&b, argv[1], "TangentKind", "Tridimensional"));
  for (it = 0; it != 7; ++it) {
    const int t = (int)types[it];
    /* mode 0: range [1, 3), mode 1: thread pool on all the points */
    for (mode = 0; mode != 2; ++mode) {
      check_status(mgis_bv_create_material_data_manager(&m, b, n));
      check_status(mgis_bv_material_data_manager_get_state_0(&s0, m));
      check_status(mgis_bv_material_data_manager_get_state_1(&s1, m));
      check_status(
          mgis_bv_material_state_manager_set_non_uniform_material_property(
              s0, "k", k, MGIS_BV_LOCAL_STORAGE));
      check_status(
          mgis_bv_material_state_manager_set_non_uniform_material_property(
              s1, "k", k, MGIS_BV_LOCAL_STORAGE));
      check_status(
          mgis_bv_material_state_manager_set_uniform_scalar_external_state_variable(
              s0, "Temperature", 293.15));
      check_status(
          mgis_bv_material_state_manager_set_uniform_scalar_external_state_variable(
              s1, "Temperature", 293.15));
      check_status(mgis_bv_material_state_manager_get_gradients(&g, s1));
      check_status(
          mgis_bv_material_state_manager_get_gradients_stride(&g_stride, s1));
      for (idx = 0; idx != n; ++idx) {
        g[idx * g_stride] = 1e-3 * (idx + 1);
      }
      if (mode == 0) {
        check_status(mgis_bv_integrate_material_data_manager_part(
            &r, m, types[it], 1, 1, 3));
      } else {
        check_status(
            mgis_bv_integrate_material_data_manager(&r, p, m, types[it], 1));
      }
      check(r == 1, "integration failed");
      check_status(mgis_bv_material_data_manager_get_tangent_operator(&K, m));
      check_status(mgis_bv_material_data_manager_get_tangent_operator_stride(
          &K_stride, m));
      check_status(mgis_bv_material_state_manager_get_internal_state_variables(
          &isvs, s1));
      check(K_stride == 36, "invalid stride of the tangent operator");
      for (idx = 0; idx != n; ++idx) {
        /* the points outside the range are left untouched */
        const int treated = (mode == 1) || ((idx >= 1) && (idx < 3));
        for (i = 0; i != 6; ++i) {
          for (j = 0; j != 6; ++j) {
            const mgis_real Kv =
                treated ? expected_operator(t, k[idx], i, j) : 0;
            check(fabs(K[idx * K_stride + i * 6 + j] - Kv) < 1e-12 * k[idx],
                  "invalid value of the operator");
          }
        }
        /* a prediction leaves the internal state variables untouched */
        check(fabs(isvs[idx] - ((treated && (t > 0)) ? 1 : 0)) < 1e-14,
              "invalid value of the internal state variable");
      }
      check_status(mgis_bv_free_material_data_manager(&m));
    }
    /* behaviour data */
    check_status(mgis_bv_allocate_behaviour_data(&d, b));
    check_status(mgis_bv_behaviour_data_get_state_0(&ds0, d));
    check_status(mgis_bv_behaviour_data_get_state_1(&ds1, d));
    check_status(mgis_bv_state_set_material_property_by_name(ds0, "k", 7));
    check_status(mgis_bv_state_set_material_property_by_name(ds1, "k", 7));
    check_status(mgis_bv_state_set_scalar_external_state_variable_by_name(
        ds0, "Temperature", 293.15));
    check_status(mgis_bv_state_set_scalar_external_state_variable_by_name(
        ds1, "Temperature", 293.15));
    check_status(mgis_bv_behaviour_data_set_time_increment(d, 1));
    check_status(mgis_bv_make_behaviour_data_view(&v, d));
    v.K[0] = t;
    check_status(mgis_bv_integrate(&r, &v, b));
    check(r == 1, "integration failed");
    check_status(mgis_bv_behaviour_data_get_tangent_operator(&K, d));
    for (i = 0; i != 6; ++i) {
      for (j = 0; j != 6; ++j) {
        check(fabs(K[i * 6 + j] - expected_operator(t, 7, i, j)) < 1e-12 * 7,
              "invalid value of the operator of the behaviour data");
      }
    }
    check_status(mgis_bv_free_behaviour_data(&d));
  }
  check_status(mgis_bv_free_behaviour(&b));
  check_status(mgis_free_thread_pool(&p));
  return test_status;
}  // end of main
