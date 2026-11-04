/*!
 * \file   tests/DerivativeTypeVariableTest.cxx
 * \brief  checks the sizes and offsets of an internal state variable whose
 * type is the derivative of a symmetric tensor type alias with respect to
 * another one
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
#include <cstdlib>
#include <iostream>
#include <string_view>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

static bool success = true;

static bool check(const bool b, const std::string_view msg) {
  if (!b) {
    success = false;
    std::cerr << msg << '\n';
  }
  return b;
}

static void check_sizes(const char* const l,
                        const mgis::behaviour::Hypothesis h) {
  using namespace mgis::behaviour;
  const auto b = load(l, "DerivativeTypeVariableTest", h);
  // derivative_type<StressStensor, StrainStensor> followed by a scalar
  const auto n = getStensorSize(h) * getStensorSize(h);
  check(b.isvs.size() == 2, "invalid number of internal state variables");
  check(getVariableSize(b.isvs.at(0), h) == n, "invalid size of 'dsig_deto'");
  check(getVariableOffset(b.isvs, "a", h) == n, "invalid offset of 'a'");
  check(getArraySize(b.isvs, h) == n + 1,
        "invalid array size for the internal state variables");
}  // end of check_sizes

int main(const int argc, const char* const* argv) {
  using namespace mgis;
  using namespace mgis::behaviour;
  constexpr auto eps = real{1e-14};
  if (!check(argc == 2, "expected two arguments")) {
    return EXIT_FAILURE;
  }
  try {
    check_sizes(argv[1], Hypothesis::TRIDIMENSIONAL);
    check_sizes(argv[1], Hypothesis::PLANESTRAIN);
    if (!success) {
      // the behaviour would write past the storage of each integration point
      return EXIT_FAILURE;
    }
    const auto b =
        load(argv[1], "DerivativeTypeVariableTest", Hypothesis::TRIDIMENSIONAL);
    auto m = MaterialDataManager{b, 2};
    setExternalStateVariable(m.s0, "Temperature", 293.15);
    setExternalStateVariable(m.s1, "Temperature", 293.15);
    const auto r = integrate(
        m, IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR, 1, 0, m.n);
    check(r == 1, "integration failed");
    const auto& isvs = m.s1.internal_state_variables;
    const auto s = m.s1.internal_state_variables_stride;
    const auto o = getVariableOffset(b.isvs, "a", b.hypothesis);
    for (size_type i = 0; i != m.n; ++i) {
      check(std::abs(isvs[i * s] - 100) < eps, "invalid value of 'dsig_deto'");
      check(std::abs(isvs[i * s + o] - 42) < eps, "invalid value of 'a'");
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
