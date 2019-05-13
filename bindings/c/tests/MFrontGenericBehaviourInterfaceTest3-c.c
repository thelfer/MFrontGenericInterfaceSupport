/*!
 * \file   MFrontGenericBehaviourInterfaceTest3-c.c
 * \brief
 * \author Thomas Helfer
 * \date   13/04/2019
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MGIS/Behaviour/Behaviour.h"

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

static void test1(const char* const* argv){
  mgis_bv_Behaviour* b;
  check_status(mgis_bv_load_behaviour(&b, argv[1],
				      argv[2], "Tridimensional"));
  mgis_size_type gs_size;
  mgis_size_type ths_size;
  check_status(
      mgis_bv_behaviour_get_gradients_size(&gs_size,b));
  check(gs_size==9,"invalid gradient size");
  check_status(
      mgis_bv_behaviour_get_thermodynamic_forces_size(&ths_size,b));
  check(ths_size==6,"invalid thermodynamic force size");
  check_status(mgis_bv_free_behaviour(&b));
}

static void test2(const char* const* argv){
  mgis_bv_Behaviour* b;
  mgis_bv_FiniteStrainBehaviourOptions* o;
  check_status(mgis_bv_create_finite_strain_behaviour_options(&o));
  check_status(mgis_bv_load_finite_strain_behaviour(&b, o, argv[1],
						    argv[2], "Tridimensional"));
  mgis_size_type gs_size;
  mgis_size_type ths_size;
  check_status(
      mgis_bv_behaviour_get_gradients_size(&gs_size,b));
  check(gs_size==9,"invalid gradient size");
  check_status(
      mgis_bv_behaviour_get_thermodynamic_forces_size(&ths_size,b));
  check(ths_size==6,"invalid thermodynamic force size");
  check_status(mgis_bv_free_finite_strain_behaviour_options(&o));
  check_status(mgis_bv_free_behaviour(&b));
}

static void test3(const char* const* argv,
		  const char* const ss, const mgis_size_type es){
  mgis_bv_Behaviour* b;
  mgis_bv_FiniteStrainBehaviourOptions* o;
  check_status(mgis_bv_create_finite_strain_behaviour_options(&o));
  check_status(mgis_bv_finite_strain_behaviour_options_set_stress_measure_by_string(o,ss));
  check_status(mgis_bv_load_finite_strain_behaviour(&b, o, argv[1],
						    argv[2], "Tridimensional"));
  mgis_size_type gs_size;
  mgis_size_type ths_size;
  check_status(
      mgis_bv_behaviour_get_gradients_size(&gs_size,b));
  check(gs_size==9,"invalid gradient size");
  check_status(
      mgis_bv_behaviour_get_thermodynamic_forces_size(&ths_size,b));
  check(ths_size==es,"invalid thermodynamic force size");
  check_status(mgis_bv_free_finite_strain_behaviour_options(&o));
  check_status(mgis_bv_free_behaviour(&b));
}

int main(const int argc, const char* const* argv) {
  if (!check(argc == 3, "expected three arguments")) {
    return EXIT_FAILURE;
  }
  test1(argv);
  test2(argv);
  test3(argv, "CAUCHY", 6);
  test3(argv, "CauchyStress", 6);
  test3(argv, "PK1", 9);
  test3(argv, "FirstPiolaKirchhoffStress", 9);
  test3(argv, "PK2", 6);
  test3(argv, "SecondPiolaKirchhoffStress", 6);
  return test_status;
}  // end of main
