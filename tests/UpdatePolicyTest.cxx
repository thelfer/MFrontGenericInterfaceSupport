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
    const auto& T_bts = m.s0.external_state_variables.at("Temperature");
    if (!std::holds_alternative<MaterialStateManager::MutableFieldHolder>(
            T_bts)) {
      std::cerr << "invalid type for the temperature\n";
      return EXIT_FAILURE;
    }
    const auto& Tvalue_bts =
        std::get<MaterialStateManager::MutableFieldHolder>(T_bts).value;
    if (!std::holds_alternative<real>(Tvalue_bts)) {
      std::cerr << "invalid type for the temperature\n";
      return EXIT_FAILURE;
    }
    const auto T = std::get<real>(Tvalue_bts);
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
    const auto& T_bts = m.s0.external_state_variables.at("Temperature");
    if (!std::holds_alternative<MaterialStateManager::MutableFieldHolder>(
            T_bts)) {
      std::cerr << "invalid type for the temperature\n";
      return EXIT_FAILURE;
    }
    const auto& Tvalue_bts =
        std::get<MaterialStateManager::MutableFieldHolder>(T_bts).value;
    if (!std::holds_alternative<real>(Tvalue_bts)) {
      std::cerr << "invalid type for the temperature\n";
      return EXIT_FAILURE;
    }
    const auto T = std::get<real>(Tvalue_bts);
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
    const auto& T_ets = m.s1.external_state_variables.at("Temperature");
    if (!std::holds_alternative<MaterialStateManager::MutableFieldHolder>(
            T_ets)) {
      std::cerr << "invalid type for the temperature\n";
      return EXIT_FAILURE;
    }
    const auto& Tvalue_ets =
        std::get<MaterialStateManager::MutableFieldHolder>(T_ets).value;
    if (!std::holds_alternative<real>(Tvalue_ets)) {
      std::cerr << "invalid type for the temperature\n";
      return EXIT_FAILURE;
    }
    const auto T = std::get<real>(Tvalue_ets);
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
    const auto& T_ets = m.s1.external_state_variables.at("Temperature");
    if (!std::holds_alternative<MaterialStateManager::MutableFieldHolder>(
            T_ets)) {
      std::cerr << "invalid type for the temperature\n";
      return EXIT_FAILURE;
    }
    const auto& Tvalue_ets =
        std::get<MaterialStateManager::MutableFieldHolder>(T_ets).value;
    if (!std::holds_alternative<real>(Tvalue_ets)) {
      std::cerr << "invalid type for the temperature\n";
      return EXIT_FAILURE;
    }
    const auto T = std::get<real>(Tvalue_ets);
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
    // a temperature bound to an external storage keeps its binding when
    // updated or reverted from a uniform value
    const auto b = load(argv[1], "Norton", Hypothesis::TRIDIMENSIONAL);
    auto check = [](const MaterialStateManager& s,
                    const std::array<real, 2>& values, const real v) {
      const auto& T = std::get<MaterialStateManager::MutableFieldHolder>(
          s.external_state_variables.at("Temperature"));
      return (std::holds_alternative<std::span<real>>(T.value)) &&
             (std::get<std::span<real>>(T.value).data() == values.data()) &&
             (std::abs(values[0] - v) < 1.e-12) &&
             (std::abs(values[1] - v) < 1.e-12);
    };
    MaterialDataManager m{b, 2};
    auto T0_values = std::array<real, 2>{T0, T0};
    setExternalStateVariable(m.s0, "Temperature", T0_values,
                             MaterialStateManager::EXTERNAL_STORAGE);
    setExternalStateVariable(m.s1, "Temperature", T1);
    update(m);
    if (!check(m.s0, T0_values, T1)) {
      std::cerr << "UpdatePolicyTest: externally stored temperature "
                   "not updated\n";
      return EXIT_FAILURE;
    }
    MaterialDataManager m2{b, 2};
    auto T1_values = std::array<real, 2>{T1, T1};
    setExternalStateVariable(m2.s0, "Temperature", T0);
    setExternalStateVariable(m2.s1, "Temperature", T1_values,
                             MaterialStateManager::EXTERNAL_STORAGE);
    revert(m2);
    if (!check(m2.s1, T1_values, T0)) {
      std::cerr << "UpdatePolicyTest: externally stored temperature "
                   "not reverted\n";
      return EXIT_FAILURE;
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
