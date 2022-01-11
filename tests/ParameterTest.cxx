/*!
 * \file   ParameterTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   13/11/2018
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
#include "MGIS/Behaviour/Behaviour.hxx"

int main(const int argc, const char* const* argv) {
  using namespace mgis::behaviour;
  constexpr const auto yg = mgis::real{150e9};
  constexpr const auto nu = mgis::real{0.3};
  constexpr const auto A = std::array<mgis::real, 2>{7, 2};
  constexpr const auto eps = mgis::real{1.e-14};
  constexpr const auto h = Hypothesis::TRIDIMENSIONAL;
  bool success = true;
  auto check = [&success](const bool b, const std::string& msg) {
    if (!b) {
      success = false;
      std::cerr << msg << '\n';
    }
    return b;
  };
  if (!check(argc == 3, "expected three arguments")) {
    return EXIT_FAILURE;
  }
  try {
    const auto d = load(argv[1], argv[2], h);
    check(d.behaviour == "ParameterTest", "invalid behaviour name");
    check(d.hypothesis == h, "invalid hypothesis");
    check(d.source == "ParameterTest.mfront", "invalid source");
    check(d.tfel_version == TFEL_VERSION, "invalid TFEL version");
    if (check(d.params.size() == 6u, "invalid number of parameters")) {
      check(d.params[0] == "YoungModulus", "invalid first parameter");
      check(d.params[1] == "PoissonRatio", "invalid second parameter");
      check(d.params[2] == "ParametersArray[0]", "invalid third parameter");
      check(d.params[3] == "ParametersArray[1]", "invalid fourth parameter");
      check(d.params[4] == "minimal_time_step_scaling_factor",
            "invalid fifth parameter");
      check(d.params[5] == "maximal_time_step_scaling_factor",
            "invalid sixth parameter");
      check(std::abs(getParameterDefaultValue<double>(d, "YoungModulus") - yg) <
                eps * yg,
            "invalid 'YoungModulus' default value");
      check(std::abs(getParameterDefaultValue<double>(d, "PoissonRatio") - nu) <
                eps * nu,
            "invalid 'PoissonRatio' default value");
      check(std::abs(getParameterDefaultValue<double>(d, "ParametersArray[0]") -
                     A[0]) < eps,
            "invalid 'ParametersArray[0]' default value");
      check(std::abs(getParameterDefaultValue<double>(d, "ParametersArray[1]") -
                     A[1]) < eps,
            "invalid 'ParametersArray[1]' default value");
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
