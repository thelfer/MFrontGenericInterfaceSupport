/*!
 * \file   IntegrateTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   03/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "MGIS/Behaviour/State.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

int main(const int argc, const char* const* argv) {
  using namespace mgis;
  using namespace mgis::behaviour;
  if (argc != 2) {
    std::cerr << "IntegrateTest: invalid number of arguments\n";
    std::exit(-1);
  }
  try {
    const auto b = load(argv[1], "NortonFailure", Hypothesis::TRIDIMENSIONAL);
    auto d = BehaviourData{b};
    const auto de = 5.e-5;
    d.dt = 180;
    // initialize the states
    setExternalStateVariable(d.s1, "Temperature", 293.15);
    // copy d.s1 in d.s0
    update(d);
    d.s1.gradients[0] = de;
    // integration
    d.K[0] = 0;
    d.rdt = 1;
    auto v = make_view(d);
    const auto r = integrate_debug(v, b);
    if ((r == 0) || (r == 1)) {
      std::cerr << "integration did not fail: " << r << "\n";
      return EXIT_FAILURE;
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
