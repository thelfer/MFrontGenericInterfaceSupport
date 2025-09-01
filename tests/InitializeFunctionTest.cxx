/*!
 * \file   tests/InitializeFunctionTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   04/02/2022
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
#include "MGIS/Behaviour/BehaviourData.hxx"
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

static void check_behaviour(const mgis::behaviour::Behaviour& b,
			    const mgis::behaviour::Hypothesis h) {
  using namespace mgis::behaviour;
  check(b.behaviour == "InitializeFunctionTest", "invalid behaviour name");
  check(b.hypothesis == h, "invalid hypothesis");
  check(b.source == "InitializeFunctionTest.mfront", "invalid source");
  check(b.tfel_version == TFEL_VERSION, "invalid TFEL version");
  check(b.initialize_functions.size() == 2u,
        "invalid number of initialize functions");
  check(b.initialize_functions.find("ElasticStrainFromInitialStress") !=
            b.initialize_functions.end(),
        "'ElasticStrainFromInitialStress' is not an initialize function");
  check(b.initialize_functions.find("StressFromInitialPressure") !=
            b.initialize_functions.end(),
        "'StressFromInitialPressure'  is not an initialize function");
  const auto& ifct1 =
      b.initialize_functions.at("ElasticStrainFromInitialStress");
  check(ifct1.inputs.empty(),
        "invalid number of inputs for the 'ElasticStrainFromInitialStress' "
        "initialize function");
  const auto& ifct2 = b.initialize_functions.at("StressFromInitialPressure");
  check(ifct2.inputs.size() == 1u,
        "invalid number of inputs for the 'StressFromInitialPressure' "
        "initialize function");
  const auto& i = ifct2.inputs.at(0);
  check(i.name == "pr",
        "invalid input name for the 'ElasticStrainFromInitialStress' "
        "initialize function");
  check(i.type == Variable::SCALAR,
        "invalid input type for the 'ElasticStrainFromInitialStress' "
        "initialize function");
  check(i.type_identifier == 0,
        "invalid input type identifier for the "
        "'ElasticStrainFromInitialStress' initialize function");
}  // end of check_behaviour

static void call_initialize_function(const mgis::behaviour::Behaviour& b) {
  using namespace mgis::behaviour;
  constexpr auto pr = mgis::real{-1.2e5};
  constexpr auto eps = -10 * pr * std::numeric_limits<mgis::real>::epsilon();
  auto d = BehaviourData{b};
  // initialize the states
  setExternalStateVariable(d.s0, "Temperature", 293.15);
  setExternalStateVariable(d.s1, "Temperature", 293.15);
  //
  auto inputs =
      allocateInitializeFunctionVariables(b, "StressFromInitialPressure");
  if (!check(inputs.size() == 1u, "invalid inputs initialisation")) {
    return;
  }
  inputs[0] = pr;
  auto v = make_view(d);
  executeInitializeFunction(v, b, "StressFromInitialPressure", inputs);
  const auto& s0 = d.s0.thermodynamic_forces;
  const auto& s1 = d.s1.thermodynamic_forces;
  if (!check((s0.size() == 6u) && (s1.size() == 6u), "invalid stress size")) {
    return;
  }
  update(d);
  for (mgis::size_type i = 0; i != 3; ++i) {
    check(std::abs(s0[i] - pr) < eps,
          "invalid stress value at the beginning of the time step");
    check(std::abs(s1[i] - pr) < eps,
          "invalid stress value at the end of the time step");
  }
  for (mgis::size_type i = 3; i != 6; ++i) {
    check(std::abs(s0[i]) < eps,
          "invalid stress value at the beginning of the time step");
    check(std::abs(s1[i]) < eps,
          "invalid stress value at the end of the time step");
  }
}  // end of call_initialize_function

static void call_initialize_function2(const mgis::behaviour::Behaviour& b) {
  using namespace mgis::behaviour;
  constexpr auto E = mgis::real{200e9};
  constexpr auto nu = mgis::real{0.3};
  constexpr auto sxx = mgis::real{150e6};
  constexpr auto eps = 10 * sxx * std::numeric_limits<mgis::real>::epsilon();
  auto d = BehaviourData{b};
  // initialize the states
  auto s_values = std::vector<mgis::real>{sxx, 0, 0, 0, 0, 0};
  d.s0.thermodynamic_forces = s_values;
  setExternalStateVariable(d.s0, "Temperature", 293.15);
  setExternalStateVariable(d.s1, "Temperature", 293.15);
  //
  auto v = make_view(d);
  executeInitializeFunction(v, b, "ElasticStrainFromInitialStress");
  update(d);
  auto eel_values =
      std::vector<mgis::real>{sxx / E, -nu * sxx / E, -nu * sxx / E, 0, 0, 0};
  const auto& s0 = d.s0.thermodynamic_forces;
  const auto& s1 = d.s1.thermodynamic_forces;
  if (!check((s0.size() == 6u) && (s1.size() == 6u), "invalid stress size")) {
    return;
  }
  const auto& eel0 = d.s0.internal_state_variables;
  const auto& eel1 = d.s1.internal_state_variables;
  for (mgis::size_type i = 0; i != 6; ++i) {
    check(std::abs(eel0[i] - eel_values[i]) < eps,
          "invalid elastic strain value at the beginning of the time step");
    check(std::abs(eel1[i] - eel_values[i]) < eps,
          "invalid elastic strain value at the beginning of the time step");
    check(std::abs(s0[i] - s_values[i]) < eps,
          "invalid stress value at the beginning of the time step");
    check(std::abs(s1[i] - s_values[i]) < eps,
          "invalid stress value at the end of the time step");
  }
}  // end of call_initialize_function

int main(const int argc, const char* const* argv) {
  using namespace mgis::behaviour;
  constexpr const auto h = Hypothesis::TRIDIMENSIONAL;
  if (!check(argc == 3, "expected three arguments")) {
    return EXIT_FAILURE;
  }
  try {
    const auto b = load(argv[1], argv[2], h);
    check_behaviour(b, h);
    call_initialize_function(b);
    call_initialize_function2(b);
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
