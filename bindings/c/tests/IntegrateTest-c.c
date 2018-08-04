/*!
 * \file   IntegrateTest-c.c
 * \brief    
 * \author Thomas Helfer
 * \date   02/08/2018
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
#include "MGIS/Behaviour/State.h"
#include "MGIS/Behaviour/BehaviourData.h"
#include "MGIS/Behaviour/BehaviourDataView.h"
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

int main(const int argc, const char* const* argv) {
  mgis_bv_State *s0, *s1;
  mgis_bv_Behaviour* b;
  mgis_bv_BehaviourData* d;
  mgis_bv_BehaviourDataView v;
  mgis_size_type o = 6;
  mgis_size_type i;
  mgis_real p[21];
  const mgis_real p_ref[21] = {0,
                               1.3523277308229e-11,
                               1.0955374667213e-07,
                               5.5890770166084e-06,
                               3.2392193670428e-05,
                               6.645865307584e-05,
                               9.9676622883138e-05,
                               0.00013302758358953,
                               0.00016635821069889,
                               0.00019969195920296,
                               0.00023302522883648,
                               0.00026635857194317,
                               0.000299691903777,
                               0.0003330252373404,
                               0.00036635857063843,
                               0.00039969190397718,
                               0.00043302523730968,
                               0.00046635857064314,
                               0.00049969190397646,
                               0.00053302523730979,
                               0.00056635857064313};
  const mgis_real de = 5.e-5;
  const mgis_real *p0, *p1;
  if (check(argc == 2, "expected two arguments") == 0) {
    return EXIT_FAILURE;
  }
  check_status(mgis_bv_load_behaviour(&b, argv[1], "Norton", "Tridimensional"));
  check_status(mgis_bv_allocate_behaviour_data(&d, b));
  // getting the offset of the equivalent plastic strain
  check_status(mgis_bv_behaviour_get_internal_state_variable_offset(
      &o, b, "EquivalentViscoplasticStrain"));
  // setting the time increment
  check_status(mgis_bv_behaviour_data_set_time_increment(d, 180));
  // state at the beginning of the time step
  check_status(mgis_bv_behaviour_data_get_state_0(&s0, d));
  // state at the end of the time step
  check_status(mgis_bv_behaviour_data_get_state_0(&s1, d));
  // setting the temperature
  check_status(mgis_bv_state_set_external_state_variable_by_name(
      s0, b, "Temperature", 293.15));
  check_status(mgis_bv_state_set_external_state_variable_by_name(
      s1, b, "Temperature", 293.15));
  // initializing the view
  check_status(mgis_bv_make_behaviour_data_view(&v, d));
  // getting the addresses where the equivalent plastic strain is stored
  check_status(mgis_bv_state_get_internal_state_variable_by_offset(&p0, s0, o));
  check_status(mgis_bv_state_get_internal_state_variable_by_offset(&p1, s1, o));
  // gradients at the end of the time step
  v.s1.gradients[0] = de;
  p[0] = *p0;
  for (i = 0; i != 20; ++i) {
    mgis_bv_integrate(&v, b);
    check_status(mgis_bv_update_behaviour_data(d));
    v.s1.gradients[0] += de;
    p[i + 1] = *p1;
  }
  for (i = 0; i != 21; ++i) {
    if (fabs(p[i] - p_ref[i]) > 1.e-12) {
      fprintf(stderr,
              "IntegrateTest-c: invalid value for "
              "the equivalent viscoplastic strain "
              "(expected '%g', computed '%g')\n",
              p_ref[i], p[i]);
      return EXIT_FAILURE;
    }
  }
  // clean-up
  check_status(mgis_bv_free_behaviour_data(&d));
  check_status(mgis_bv_free_behaviour(&b));
  return test_status;
} // end of main
