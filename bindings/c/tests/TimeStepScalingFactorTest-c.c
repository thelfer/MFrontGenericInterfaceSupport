/*!
 * \file   TimeStepScalingFactorTest-c.c
 * \brief  tests of the time step scaling factor returned by the integrations
 * of a material data manager
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
#include "MGIS/Behaviour/MaterialDataManager.h"
#include "MGIS/Behaviour/Integrate.h"

int test_status = EXIT_SUCCESS;
mgis_bv_Behaviour* b = NULL;
mgis_bv_MaterialDataManager* m = NULL;
mgis_bv_BehaviourIntegrationOptions* o = NULL;
mgis_ThreadPool* p = NULL;

static void free_all(void) {
  mgis_bv_free_behaviour_integration_options(&o);
  mgis_bv_free_material_data_manager(&m);
  mgis_bv_free_behaviour(&b);
  mgis_free_thread_pool(&p);
}  // end of free_all

static void check_status(const mgis_status s) {
  if (s.exit_status != MGIS_SUCCESS) {
    fprintf(stderr, "invalid function call: %s\n", s.msg);
    free_all();
    exit(EXIT_FAILURE);
  }
}  // end of check_status

static int check(const int c, const char* const e) {
  if (c == 0) {
    test_status = EXIT_FAILURE;
    fprintf(stderr, "%s\n", e);
  }
  return c;
}  // end of check

int main(const int argc, const char* const* argv) {
  /* time step scaling factors proposed at the three integration points */
  mgis_real factors[3] = {2, 0.5, 3};
  mgis_bv_MaterialStateManager* s1;
  mgis_real tsf;
  int r;
  if (check(argc == 2, "expected two arguments") == 0) {
    return EXIT_FAILURE;
  }
  check_status(mgis_create_thread_pool(&p, 2));
  check_status(mgis_bv_load_behaviour(&b, argv[1], "TimeStepScalingFactorTest",
                                      "Tridimensional"));
  check_status(mgis_bv_create_material_data_manager(&m, b, 3));
  check_status(mgis_bv_material_data_manager_get_state_1(&s1, m));
  check_status(mgis_bv_material_state_manager_set_non_uniform_material_property(
      s1, "ProposedTimeStepScalingFactor", factors, MGIS_BV_LOCAL_STORAGE));
  check_status(
      mgis_bv_material_state_manager_set_uniform_scalar_external_state_variable(
          s1, "Temperature", 293.15));
  check_status(mgis_bv_update_material_data_manager(m));
  check_status(mgis_bv_create_behaviour_integration_options(&o));
  check_status(mgis_bv_behaviour_integration_options_set_integration_type(
      o, MGIS_BV_INTEGRATION_NO_TANGENT_OPERATOR));
  /* all the integration points: minimum of the proposed factors */
  check_status(mgis_bv_integrate_material_data_manager_with_options_2(
      &r, &tsf, p, m, o, 1));
  check(fabs(tsf - 0.5) < 1e-14, "invalid time step scaling factor");
  check_status(mgis_bv_integrate_debug_material_data_manager_with_options_2(
      &r, &tsf, p, m, o, 1));
  check(fabs(tsf - 0.5) < 1e-14, "invalid time step scaling factor");
  /* last integration point only */
  check_status(mgis_bv_integrate_material_data_manager_part_2(
      &r, &tsf, m, MGIS_BV_INTEGRATION_NO_TANGENT_OPERATOR, 1, 2, 3));
  check(r == 1, "invalid result");
  check(fabs(tsf - 3) < 1e-14, "invalid time step scaling factor");
  check_status(mgis_bv_integrate_debug_material_data_manager_part_2(
      &r, &tsf, m, MGIS_BV_INTEGRATION_NO_TANGENT_OPERATOR, 1, 2, 3));
  check(fabs(tsf - 3) < 1e-14, "invalid time step scaling factor");
  free_all();
  return test_status;
}  // end of main
