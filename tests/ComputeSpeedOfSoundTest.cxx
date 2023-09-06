/*!
 * \file   tests/ComputeSpeedOfSoundTest.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   06/09/2023
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
    const auto b = load(argv[1], "ComputeSpeedOfSound", Hypothesis::TRIDIMENSIONAL);
    auto d = BehaviourData{b};
    const auto de = 5.e-5;
    d.dt = 180;
    // initialize the states
    setExternalStateVariable(d.s1, "Temperature", 293.15);
    d.s1.mass_density = 7850;
    // copy d.s1 in d.s0
    update(d);
    d.s1.gradients[0] = de;
    // integration
    d.K[0] = 100;
    d.rdt = 1;
    auto v = make_view(d);
    integrate(v, b);
    update(d);
    const auto yg = 150e9;
    const auto v_ref = sqrt(yg / d.s1.mass_density);
    if (std::abs(d.speed_of_sound - v_ref) > 1.e-12 * v_ref) {
      std::cerr << "IntegrateTest: invalid value for the speed of sound"
                << "(expected '" << v_ref << "', computed '" << d.speed_of_sound
                << "')\n";
      return EXIT_FAILURE;
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}


