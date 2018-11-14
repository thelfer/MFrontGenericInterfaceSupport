/*!
 * \file   BoundsCheckTest.cxx
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
  using lreal = long double;
  constexpr const auto eps = lreal{1.e-14};
  constexpr const auto yg_min = lreal{100.e9};
  constexpr const auto yg_max = lreal{200.e9};
  constexpr const auto nu_min = lreal{-1};
  constexpr const auto nu_max = lreal{0.5};
  constexpr const auto iv_pmin = lreal{0};
  constexpr const auto iv_pmax = lreal{0.8};
  constexpr const auto iv_min = lreal{0.2};
  constexpr const auto iv_max = lreal{0.5};
  constexpr const auto ev_pmin = lreal{0};
  constexpr const auto ev_pmax = lreal{500};
  constexpr const auto ev_min = lreal{200};
  constexpr const auto ev_max = lreal{400};
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
    check(d.behaviour == "BoundsCheckTest", "invalid behaviour name");
    check(d.hypothesis == h, "invalid hypothesis");
    check(d.source == "BoundsCheckTest.mfront", "invalid source");
    check(d.tfel_version == TFEL_VERSION, "invalid TFEL version");
    // test on material properties
    check(hasBounds(d, "YoungModulus"), "'YoungModulus' shall have bounds");
    check(hasLowerBound(d, "YoungModulus"),
          "'YoungModulus' shall have a lower bound");
    check(std::abs(getLowerBound(d, "YoungModulus") - yg_min) < eps * yg_min,
          "invalid value for the 'YoungModulus' lower bound");
    check(std::abs(getUpperBound(d, "YoungModulus") - yg_max) < eps * yg_max,
          "invalid value for the 'YoungModulus' upper bound");
    check(hasUpperBound(d, "YoungModulus"),
          "'YoungModulus' shall have an upper bound");
    check(hasPhysicalBounds(d, "YoungModulus"),
          "'YoungModulus' shall have physical bounds");
    check(hasLowerPhysicalBound(d, "YoungModulus"),
          "'YoungModulus' shall have a lower physical bound");
    check(std::abs(getLowerPhysicalBound(d, "YoungModulus")) < eps * yg_max,
          "invalid value for the 'YoungModulus' physical lower bound");
    check(!hasUpperPhysicalBound(d, "YoungModulus"),
          "'YoungModulus' shall not have an upper physical bound");
    check(!hasBounds(d, "PoissonRatio"),
          "'PoissonRatio' shall not have bounds");
    check(!hasLowerBound(d, "PoissonRatio"),
          "'PoissonRatio' shall not have a lower bound");
    check(!hasUpperBound(d, "PoissonRatio"),
          "'PoissonRatio' shall not have an upper bound");
    check(hasPhysicalBounds(d, "PoissonRatio"),
          "'PoissonRatio' shall have physical bounds");
    check(hasLowerPhysicalBound(d, "PoissonRatio"),
          "'PoissonRatio' shall have a lower physical bound");
    check(std::abs(getLowerPhysicalBound(d, "PoissonRatio") - nu_min) < eps,
          "invalid value for the 'PoissonRatio' physical lower bound");
    check(hasUpperPhysicalBound(d, "PoissonRatio"),
          "'PoissonRatio' shall have an upper physical bound");
    check(std::abs(getUpperPhysicalBound(d, "PoissonRatio") - nu_max) <
              eps * nu_max,
          "invalid value for the 'PoissonRatio' physical upper bound");
    // internal state variables
    check(hasBounds(d, "StateVariable"), "'StateVariable' shall have bounds");
    check(hasLowerBound(d, "StateVariable"),
          "'StateVariable' shall have a lower bound");
    check(std::abs(getLowerBound(d, "StateVariable") - iv_min) < eps * iv_min,
          "invalid value for the 'StateVariable' lower bound");
    check(std::abs(getUpperBound(d, "StateVariable") - iv_max) < eps * iv_max,
          "invalid value for the 'StateVariable' upper bound");
    check(hasUpperBound(d, "StateVariable"),
          "'StateVariable' shall have an upper bound");
    check(hasPhysicalBounds(d, "StateVariable"),
          "'StateVariable' shall have physical bounds");
    check(hasLowerPhysicalBound(d, "StateVariable"),
          "'StateVariable' shall have a lower physical bound");
    check(std::abs(getLowerPhysicalBound(d, "StateVariable") - iv_pmin) < eps,
          "invalid value for the 'StateVariable' physical lower bound");
    check(hasUpperPhysicalBound(d, "StateVariable"),
          "'StateVariable' shall not have an upper physical bound");
    check(std::abs(getUpperPhysicalBound(d, "StateVariable") - iv_pmax) <
              eps * iv_pmax,
          "invalid value for the 'StateVariable' physical upper bound");
    // external state variable
    check(hasBounds(d, "ExternalStateVariable"),
          "'ExternalStateVariable' shall have bounds");
    check(hasLowerBound(d, "ExternalStateVariable"),
          "'ExternalStateVariable' shall have a lower bound");
    check(std::abs(getLowerBound(d, "ExternalStateVariable") - ev_min) <
              eps * ev_min,
          "invalid value for the 'ExternalStateVariable' lower bound");
    check(std::abs(getUpperBound(d, "ExternalStateVariable") - ev_max) <
              eps * ev_max,
          "invalid value for the 'ExternalStateVariable' upper bound");
    check(hasUpperBound(d, "ExternalStateVariable"),
          "'ExternalStateVariable' shall have an upper bound");
    check(hasPhysicalBounds(d, "ExternalStateVariable"),
          "'ExternalStateVariable' shall have physical bounds");
    check(hasLowerPhysicalBound(d, "ExternalStateVariable"),
          "'ExternalStateVariable' shall have a lower physical bound");
    check(std::abs(getLowerPhysicalBound(d, "ExternalStateVariable") -
                   ev_pmin) < eps,
          "invalid value for the 'ExternalStateVariable' physical lower bound");
    check(hasUpperPhysicalBound(d, "ExternalStateVariable"),
          "'ExternalStateVariable' shall not have an upper physical bound");
    check(std::abs(getUpperPhysicalBound(d, "ExternalStateVariable") -
                   ev_pmax) < eps * ev_pmax,
          "invalid value for the 'ExternalStateVariable' physical upper bound");
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
