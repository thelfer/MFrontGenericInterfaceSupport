/*!
 * \file   tests/FiniteStrainPostProcessingTest.cxx
 * \brief  checks that the post-processings of finite strain behaviours use
 * the stress measure chosen at loading
 * \date   20/10/2026
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <cmath>
#include <limits>
#include <vector>
#include <cstdlib>
#include <iostream>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

int main(const int argc, const char* const* argv) {
  using namespace mgis;
  using namespace mgis::behaviour;
  if (argc != 2) {
    std::cerr
        << "FiniteStrainPostProcessingTest: invalid number of arguments\n";
    return EXIT_FAILURE;
  }
  try {
    auto o = FiniteStrainBehaviourOptions{};
    o.stress_measure = FiniteStrainBehaviourOptions::PK1;
    o.tangent_operator = FiniteStrainBehaviourOptions::DPK1_DF;
    const auto b = load(o, argv[1], "FiniteStrainPostProcessingTest",
                        Hypothesis::TRIDIMENSIONAL);
    auto m = MaterialDataManager{b, 1};
    setExternalStateVariable(m.s0, "Temperature", 293.15);
    setExternalStateVariable(m.s1, "Temperature", 293.15);
    const real F[9] = {1.01, 1, 1, 0, 0, 0, 0, 0, 0};
    std::copy(F, F + 9, m.s1.gradients.begin());
    const auto r = integrate(
        m, IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR, 1, 0, m.n);
    if (r != 1) {
      std::cerr << "integration failed\n";
      return EXIT_FAILURE;
    }
    auto pr = std::vector<real>(1);
    const auto rp = executePostProcessing(pr, m, "MeanStress");
    if (rp.exit_status != 1) {
      std::cerr << "post-processing failed\n";
      return EXIT_FAILURE;
    }
    // mean of the Kirchhoff stress, computed from the first Piola-Kirchhoff
    // stress, the deformation gradient being diagonal. The Kirchhoff stress
    // is equal to the dual of the Hencky strain for this isotropic behaviour.
    const auto& P = m.s1.thermodynamic_forces;
    const auto pr_ref = (P[0] * F[0] + P[1] * F[1] + P[2] * F[2]) / 3;
    if (std::abs(pr[0] - pr_ref) >
        100 * std::numeric_limits<real>::epsilon() * std::abs(pr_ref)) {
      std::cerr << "invalid mean stress: " << pr[0] << " vs " << pr_ref << '\n';
      return EXIT_FAILURE;
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
