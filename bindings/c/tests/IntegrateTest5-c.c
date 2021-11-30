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
#include "MGIS/Behaviour/BehaviourData.h"
#include "MGIS/Behaviour/Integrate.h"

int test_status = EXIT_SUCCESS;
mgis_model_Model* model = NULL;
mgis_bv_BehaviourData* d = NULL;

static void check_status(const mgis_status s) {
  if (s.exit_status != MGIS_SUCCESS) {
    fprintf(stderr, "invalid function call: %s\n", s.msg);
    mgis_model_free_model(&model);
    mgis_bv_free_behaviour_data(&d);
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
  mgis_bv_BehaviourDataView v;
  mgis_size_type o;
  mgis_real* isvs0; /* internal state variables at the beginning of the time step */
  mgis_real* isvs1; /* internal state variables at the end of the time step */
  mgis_bv_State* s0; /* state at the beginning of the time step */
  mgis_bv_State* s1; /* state at the end of the time step */
  mgis_real xvalues[11]; /* values of 'x' for the first integration point */
  mgis_real* x;
  mgis_size_type idx;
  mgis_size_type i;
  mgis_real A;
  mgis_real t;
  mgis_real x_ref;
  int r;
  check_status(mgis_model_load(&model, argv[1], "ode_rk54", "Tridimensional"));
  check_status(mgis_bv_behaviour_get_parameter_default_value(&A, model, "A"));
  check_status(mgis_bv_allocate_behaviour_data(&d, model));
  check_status(mgis_bv_behaviour_get_internal_state_variable_offset_by_name(
      &o, model, "x"));
  check_status(mgis_bv_behaviour_data_get_state_0(&s0, d));
  check_status(mgis_bv_behaviour_data_get_state_1(&s1, d));
  /* initialize the internal state variable */
  mgis_bv_state_get_internal_state_variable_by_offset(&x, s1, o);
  *x = 1;
  /* initialize the external state variable */
  check_status(mgis_bv_state_set_scalar_external_state_variable_by_name(
      s1, "Temperature", 293.15));
  /* copy s1 in s0 */
  check_status(mgis_bv_update_behaviour_data(d));
  // setting the time increment
  check_status(mgis_bv_behaviour_data_set_time_increment(d, dt));
  // initializing the view
  check_status(mgis_bv_make_behaviour_data_view(&v, d));
  // integration */
  xvalues[0] = *x;
  for (i = 0; i != 10; ++i) {
    v.K[0] = 0;
    check_status(mgis_bv_integrate(&r, &v, model));
    check_status(mgis_bv_update_behaviour_data(d));
    xvalues[i+1] = *x;
  }
  // clean-up
  check_status(mgis_model_free_model(&model));
  check_status(mgis_bv_free_behaviour_data(&d));
  // checks
  t = 0;
  for (i = 0; i != 11; ++i) {
    x_ref = exp(-A * t);
    if (fabs(xvalues[i] - x_ref) > eps) {
      fprintf(stderr,
              "IntegrateTest: invalid value for x "
              "at the first integration point"
              "(expected '%g', computed '%g')\n",
              x_ref, xvalues[i]);
      return EXIT_FAILURE;
    } 
    t += dt;
  }
  return EXIT_SUCCESS;
}
