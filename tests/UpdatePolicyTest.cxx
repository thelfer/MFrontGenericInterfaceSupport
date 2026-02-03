/*!
 * \file   UpdatePolicyTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   25/01/2026
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
  constexpr auto T0 = real{293.15};
  constexpr auto T1 = real{893.15};
  if (argc != 2) {
    std::cerr << "IntegrateTest: invalid number of arguments\n";
    std::exit(-1);
  }
  try {
    const auto b = load(argv[1], "Norton", Hypothesis::TRIDIMENSIONAL);
    MaterialDataManager m{b, 100};
    // initialize the external state variable
    setExternalStateVariable(m.s0, "Temperature", T0,
                             MaterialStateManager::UPDATE);
    setExternalStateVariable(m.s1, "Temperature", T1,
                             MaterialStateManager::UPDATE);
    update(m);
    const auto T =
        std::get<real>(m.s0.external_state_variables["Temperature"].value);
    if (std::abs(T - T1) > 1.e-12) {
      std::cerr << "UpdatePolicyTest: invalid temperature value (" << T
                << ")\n";
      return EXIT_FAILURE;
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  try {
    const auto b = load(argv[1], "Norton", Hypothesis::TRIDIMENSIONAL);
    MaterialDataManager m{b, 100};
    // initialize the external state variable
    setExternalStateVariable(m.s0, "Temperature", T0,
                             MaterialStateManager::NOUPDATE);
    setExternalStateVariable(m.s1, "Temperature", T1,
                             MaterialStateManager::NOUPDATE);
    update(m);
    const auto T =
        std::get<real>(m.s0.external_state_variables["Temperature"].value);
    if (std::abs(T - T0) > 1.e-12) {
      std::cerr << "UpdatePolicyTest: invalid temperature value (" << T
                << ")\n";
      return EXIT_FAILURE;
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  try {
    const auto b = load(argv[1], "Norton", Hypothesis::TRIDIMENSIONAL);
    MaterialDataManager m{b, 100};
    // initialize the external state variable
    setExternalStateVariable(m.s0, "Temperature", T0,
                             MaterialStateManager::UPDATE);
    setExternalStateVariable(m.s1, "Temperature", T1,
                             MaterialStateManager::UPDATE);
    revert(m);
    const auto T =
        std::get<real>(m.s1.external_state_variables["Temperature"].value);
    if (std::abs(T - T0) > 1.e-12) {
      std::cerr << "UpdatePolicyTest: invalid temperature value (" << T
                << ")\n";
      return EXIT_FAILURE;
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  try {
    const auto b = load(argv[1], "Norton", Hypothesis::TRIDIMENSIONAL);
    MaterialDataManager m{b, 100};
    // initialize the external state variable
    setExternalStateVariable(m.s0, "Temperature", T0,
                             MaterialStateManager::NOUPDATE);
    setExternalStateVariable(m.s1, "Temperature", T1,
                             MaterialStateManager::NOUPDATE);
    revert(m);
    const auto T =
        std::get<real>(m.s1.external_state_variables["Temperature"].value);
    if (std::abs(T - T1) > 1.e-12) {
      std::cerr << "UpdatePolicyTest: invalid temperature value (" << T
                << ")\n";
      return EXIT_FAILURE;
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
