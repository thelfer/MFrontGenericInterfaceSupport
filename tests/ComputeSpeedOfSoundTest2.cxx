/*!
 * \file   ComputeSpeedOfSoundTest2.cxx
 * \brief
 * \author Thomas Helfer
 * \date   24/08/2018
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
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

int main(const int argc, const char* const* argv) {
  using namespace mgis;
  using namespace mgis::behaviour;
  const auto yg = 150e9;
  const auto rho = real{7850};
  if (argc != 2) {
    std::cerr << "ComputeSpeedOfSoundTest: invalid number of arguments\n";
    std::exit(-1);
  }
  try {
    const auto b = load(argv[1], "ComputeSpeedOfSound", Hypothesis::TRIDIMENSIONAL);
    MaterialDataManager m{b, 100};
    const auto de = 5.e-5;
    // initialize the external state variable
    m.s1.external_state_variables["Temperature"] = 293.15;
    m.s1.mass_density = rho;
    // copy d.s1 in d.s0
    update(m);
    for (size_type idx = 0; idx != m.n; ++idx) {
      m.s1.gradients[idx * m.s1.gradients_stride] = de;
    }
    const auto dt = real(180);
    auto opts = BehaviourIntegrationOptions{};
    opts.integration_type = IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR;
    opts.compute_speed_of_sound = true;
    integrate(m, opts, dt, 0, m.n);
    std::cerr.precision(14);
    const auto v_ref = std::sqrt(yg/rho);
    if (std::abs(m.speed_of_sound[0] - v_ref) > 1.e-12 * v_ref) {
      std::cerr
          << "ComputeSpeedOfSoundTest: invalid value for the speed of sound\n";
      return EXIT_FAILURE;
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
