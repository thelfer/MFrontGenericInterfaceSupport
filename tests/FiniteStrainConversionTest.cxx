/*!
 * \file   tests/FiniteStrainConversionTest.cxx
 * \brief  tests of the conversions of the Cauchy stress and of its
 * derivative with respect to the deformation gradient to the first
 * Piola-Kirchhoff stress and to its derivative: the converted values are
 * compared to the ones returned by the same behaviour loaded with the
 * `PK1` and `DPK1_DF` options.
 * \date   03/11/2026
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <cmath>
#include <span>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/Integrate.hxx"
#include "MGIS/Behaviour/FiniteStrainSupport.hxx"

static bool success = true;

static void check(const bool b, const char* const msg) {
  if (!b) {
    std::cerr << "FiniteStrainConversionTest: " << msg << '\n';
    success = false;
  }
}  // end of check

//! \return true if the given arrays are equal up to a relative tolerance
static bool are_equal(std::span<const mgis::real> a,
                      std::span<const mgis::real> b) {
  if (a.size() != b.size()) {
    return false;
  }
  auto m = mgis::real{0};
  auto e = mgis::real{0};
  for (std::size_t i = 0; i != b.size(); ++i) {
    m = std::max(m, std::abs(b[i]));
    e = std::max(e, std::abs(a[i] - b[i]));
  }
  return e < 1.e-10 * m;
}  // end of are_equal

int main(const int argc, const char* const* argv) {
  using namespace mgis;
  using namespace mgis::behaviour;
  if (argc != 2) {
    std::cerr << "FiniteStrainConversionTest: invalid number of arguments\n";
    std::exit(-1);
  }
  try {
    auto oc = FiniteStrainBehaviourOptions{};
    oc.stress_measure = FiniteStrainBehaviourOptions::CAUCHY;
    oc.tangent_operator = FiniteStrainBehaviourOptions::DSIG_DF;
    auto op = FiniteStrainBehaviourOptions{};
    op.stress_measure = FiniteStrainBehaviourOptions::PK1;
    op.tangent_operator = FiniteStrainBehaviourOptions::DPK1_DF;
    for (const auto h :
         {Hypothesis::TRIDIMENSIONAL, Hypothesis::AXISYMMETRICAL}) {
      const auto bc = load(oc, argv[1], "FiniteStrainPostProcessingTest", h);
      const auto bp = load(op, argv[1], "FiniteStrainPostProcessingTest", h);
      const auto ts = getTensorSize(h);
      // deformation gradients with shear and a non-unit determinant,
      // different at each point
      const auto n = size_type{3};
      auto F = std::vector<real>(n * ts);
      for (size_type i = 0; i != n; ++i) {
        const auto a = 0.01 * (i + 1);
        const real F3D[9] = {1 + a,    1 - a / 2, 1 + 2 * a, 0.3 * a, -0.2 * a,
                             0.15 * a, 0.25 * a,  -0.1 * a,  0.4 * a};
        std::copy(F3D, F3D + ts, F.begin() + ts * i);
      }
      // material data managers
      MaterialDataManager mc(bc, n), mp(bp, n);
      for (auto* const m : {&mc, &mp}) {
        setExternalStateVariable(m->s0, "Temperature", 293.15);
        setExternalStateVariable(m->s1, "Temperature", 293.15);
        std::copy(F.begin(), F.end(), m->s1.gradients.begin());
        const auto r = integrate(
            *m, IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR, 1, 0,
            n);
        check(r == 1, "integration failed");
      }
      auto P = std::vector<real>(n * ts);
      auto dP = std::vector<real>(n * ts * ts);
      auto Ps = std::span<real>(P);
      auto dPs = std::span<real>(dP);
      convertFiniteStrainStress(Ps, mc, FiniteStrainStress::PK1);
      convertFiniteStrainTangentOperator(dPs, mc,
                                         FiniteStrainTangentOperator::DPK1_DF);
      check(are_equal(P, mp.s1.thermodynamic_forces),
            "invalid conversion of the stress (material data manager)");
      check(are_equal(dP, mp.K),
            "invalid conversion of the tangent operator "
            "(material data manager)");
      // behaviour data, deformation gradient of the last point
      BehaviourData dc(bc), dp(bp);
      for (auto* const d : {&dc, &dp}) {
        setExternalStateVariable(d->s0, "Temperature", 293.15);
        setExternalStateVariable(d->s1, "Temperature", 293.15);
        std::copy(F.begin() + (n - 1) * ts, F.end(), d->s1.gradients.begin());
        d->dt = 1;
        d->K[0] = static_cast<int>(
            IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR);
        auto v = make_view(*d);
        check(integrate(v, d == &dc ? bc : bp) == 1, "integration failed");
      }
      auto P1 = std::vector<real>(ts);
      auto dP1 = std::vector<real>(ts * ts);
      auto P1s = std::span<real>(P1);
      auto dP1s = std::span<real>(dP1);
      convertFiniteStrainStress(P1s, dc, FiniteStrainStress::PK1);
      convertFiniteStrainTangentOperator(dP1s, dc,
                                         FiniteStrainTangentOperator::DPK1_DF);
      check(are_equal(P1, dp.s1.thermodynamic_forces),
            "invalid conversion of the stress (behaviour data)");
      check(are_equal(dP1, dp.K),
            "invalid conversion of the tangent operator (behaviour data)");
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
