/*!
 * \file   IntegrateTest2-c.c
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
#include "MGIS/Behaviour/MaterialDataManager.h"
#include "MGIS/Behaviour/Integrate.h"

int test_status = EXIT_SUCCESS;
mgis_bv_Behaviour* b = NULL;
mgis_bv_MaterialDataManager* m = NULL;

static void check_status(const mgis_status s) {
  if (s.exit_status != MGIS_SUCCESS) {
    fprintf(stderr, "invalid function call: %s\n", s.msg);
    mgis_bv_free_behaviour(&b);
    mgis_bv_free_material_data_manager(&m);
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
  mgis_size_type o;
  if (check(argc == 2, "expected two arguments") == 0) {
    return EXIT_FAILURE;
  }
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
  const mgis_real de = 5.e-5; /* strain increment */
  const mgis_real dt = 180;
  mgis_real* g; /* gradients */
  mgis_size_type g_stride; /* gradients size */
  mgis_real* isvs0; /* internal state variables at the beginning of the time step */
  mgis_real* isvs1; /* internal state variables at the end of the time step */
  mgis_size_type isvs_stride; /* internal state variables stride */
  const mgis_size_type n = 100;
  mgis_bv_MaterialStateManager*
      s0; /* state at the beginning of the time step */
  mgis_bv_MaterialStateManager* s1; /* state at the end of the time step */
  mgis_real pi[21]; /* values of the equivalent plastic strain for the first
                       integration point */
  mgis_real pe[21]; /* values of the equivalent plastic strain for the last
                       integration point */
  mgis_size_type ni,ne;
  mgis_size_type idx;
  mgis_size_type i;
  int r;
  check_status(mgis_bv_load_behaviour(&b, argv[1], "Norton", "Tridimensional"));
  check_status(mgis_bv_create_material_data_manager(&m,b,100));
  check_status(mgis_bv_behaviour_get_internal_state_variable_offset(
      &o, b, "EquivalentViscoplasticStrain"));
  check_status(mgis_bv_behaviour_get_internal_state_variable_offset(
      &o, b, "EquivalentViscoplasticStrain"));
  check_status(mgis_bv_material_data_manager_get_state_0(&s0, m));
  check_status(mgis_bv_material_data_manager_get_state_1(&s1, m));
  /* initialize the external state variable */
  check_status(
      mgis_bv_material_state_manager_set_uniform_external_state_variable(
          s1, "Temperature", 293.15));
  /* copy s1 in s0 */
  check_status(mgis_bv_update_material_data_manager(m));
  check_status(mgis_bv_material_state_manager_get_gradients(&g,s1));
  check_status(mgis_bv_material_state_manager_get_gradients_stride(&g_stride,s1));
  check_status(
      mgis_bv_material_state_manager_get_internal_state_variables(&isvs0, s0));
  check_status(
      mgis_bv_material_state_manager_get_internal_state_variables(&isvs1, s1));
  check_status(
      mgis_bv_material_state_manager_get_internal_state_variables_stride(
          &isvs_stride, s1));
  for (idx = 0; idx != n; ++idx) {
    g[idx * g_stride] = de;
  }
  ni = o;
  ne = (n - 1) * isvs_stride + o;
  /* // integration */
  pi[0] = isvs0[ni];
  pe[0] = isvs0[ne];
  for (i = 0; i != 20; ++i) {
    check_status(mgis_bv_integrate_material_data_manager_part(&r, m, dt, 0, n));
    check_status(mgis_bv_update_material_data_manager(m));
    for (idx = 0; idx != n; ++idx) {
      g[idx * g_stride] += de;
    }
    pi[i + 1] = isvs1[ni];
    pe[i + 1] = isvs1[ne];
  }
  check_status(mgis_bv_free_behaviour(&b));
  check_status(mgis_bv_free_material_data_manager(&m));
  for (i = 0; i != 21; ++i) {
    if (fabs(pi[i] - p_ref[i]) > 1.e-12) {
      fprintf(stderr,
              "IntegrateTest: invalid value for the equivalent "
              "viscoplastic strain at the first integration point"
              "(expected '%g', computed '%g')\n",
              p_ref[i], pi[i]);
      return EXIT_FAILURE;
    } 
    if (fabs(pe[i] - p_ref[i]) > 1.e-12) {
      fprintf(stderr,
              "IntegrateTest: invalid value for the equivalent "
              "viscoplastic strain at the first integration point"
              "(expected '%g', computed '%g')\n",
              p_ref[i], pe[i]);
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}
