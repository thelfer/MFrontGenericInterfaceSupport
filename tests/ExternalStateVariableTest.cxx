/*!
 * \file   tests/ExternalStateVariableTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   29/11/2021
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <span>
#include <array>
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
  check(b.behaviour == "TensorialExternalStateVariableTest",
        "invalid behaviour name");
  check(b.hypothesis == h, "invalid hypothesis");
  check(b.source == "TensorialExternalStateVariableTest.mfront",
        "invalid source");
  check(b.tfel_version == TFEL_VERSION, "invalid TFEL version");
  check(b.esvs.size() == 10, "invalid number of external state variables");
  const auto esvnames = std::vector<std::string>{
      "Temperature", "v_esv",     "v2_esv[0]", "v2_esv[1]", "s_esv",
      "s2_esv[0]",   "s2_esv[1]", "t_esv",     "t2_esv[0]", "t2_esv[1]"};
  const auto esvtypes = std::vector<Variable::Type>{
      Variable::SCALAR,  Variable::VECTOR,  Variable::VECTOR,  Variable::VECTOR,
      Variable::STENSOR, Variable::STENSOR, Variable::STENSOR, Variable::TENSOR,
      Variable::TENSOR,  Variable::TENSOR};
  const auto esvsizes =
      std::vector<mgis::size_type>{1, 3, 3, 3, 6, 6, 6, 9, 9, 9};
  for (mgis::size_type i = 0; i != 10; ++i) {
    check(b.esvs.at(i).name == esvnames.at(i),
          "invalid external state variable name");
    check(b.esvs.at(i).type == esvtypes.at(i),
          "invalid external variable type");
    check(getVariableSize(b.esvs.at(i), b.hypothesis) == esvsizes.at(i),
          "invalid external variable type");
  }
  check(getArraySize(b.esvs, b.hypothesis) == 55,
        "invalid array size for the external state variables");
}  // end of check_behaviour

static void check_behaviour_data(const mgis::behaviour::Behaviour& b) {
  using namespace mgis::behaviour;
  constexpr const auto eps = 1.e-14;
  BehaviourData d{b};
  check(d.s0.external_state_variables.size() == 55,
        "invalid array size for the external state variables");
  check(d.s1.external_state_variables.size() == 55,
        "invalid array size for the external state variables");
  const auto v_esv_values = std::vector<mgis::real>{1, 2, 3};
  setExternalStateVariable(d.s1, "v_esv", v_esv_values);
  auto view = make_view(d);
  integrate(view, b);
  for (mgis::size_type i = 0; i != 3; ++i) {
    check(std::abs(d.s1.external_state_variables.at(i + 1) -
                   v_esv_values.at(i)) < eps,
          "invalid external state variable value");
    check(std::abs(d.s1.internal_state_variables.at(i) - v_esv_values.at(i)) <
              eps,
          "invalid internal state variable value");
  }
}  // end of check_behaviour_data

static void check_material_data_manager(const mgis::behaviour::Behaviour& b) {
  using namespace mgis::behaviour;
  constexpr const auto eps = 1.e-14;
  constexpr const auto dt = 0.1;
  MaterialDataManager d{b, 2};
  setMaterialProperty(d.s0, "YoungModulus", 150e9);
  setMaterialProperty(d.s0, "PoissonRatio", 0.3);
  setMaterialProperty(d.s1, "YoungModulus", 150e9);
  setMaterialProperty(d.s1, "PoissonRatio", 0.3);
  auto zeros = std::vector<mgis::real>(9, 0);
  auto v_esv_values = std::vector<mgis::real>{1, 2, 3};
  // v_esv is uniform
  setExternalStateVariable(d.s0, "v_esv",
                           std::span<mgis::real>(zeros.data(), 3u),
                           MaterialStateManager::EXTERNAL_STORAGE);
  setExternalStateVariable(d.s1, "v_esv", v_esv_values);
  // v2_esv[0] is not uniform
  auto v2_esv0_values = std::vector<mgis::real>{1, 2, 3, 7, 6, 5};
  setExternalStateVariable(d.s0, "v2_esv[0]",
                           std::span<mgis::real>(zeros.data(), 3u),
                           MaterialStateManager::EXTERNAL_STORAGE);
  setExternalStateVariable(d.s1, "v2_esv[0]", v2_esv0_values);
  auto v2_esv1_values = std::vector<mgis::real>{9, 8, 2, 3, 6, 4};
  setExternalStateVariable(d.s0, "v2_esv[1]",
                           std::span<mgis::real>(zeros.data(), 3u),
                           MaterialStateManager::EXTERNAL_STORAGE);
  setExternalStateVariable(d.s1, "v2_esv[1]", v2_esv1_values,
                           MaterialStateManager::EXTERNAL_STORAGE);
  auto init = [&zeros](MaterialStateManager& s) {
    setExternalStateVariable(s, "Temperature",
                             std::span<mgis::real>(zeros.data(), 1u),
                             MaterialStateManager::EXTERNAL_STORAGE);
    setExternalStateVariable(s, "s_esv",
                             std::span<mgis::real>(zeros.data(), 6u),
                             MaterialStateManager::EXTERNAL_STORAGE);
    setExternalStateVariable(s, "s2_esv[0]",
                             std::span<mgis::real>(zeros.data(), 6u),
                             MaterialStateManager::EXTERNAL_STORAGE);
    setExternalStateVariable(s, "s2_esv[1]",
                             std::span<mgis::real>(zeros.data(), 6u),
                             MaterialStateManager::EXTERNAL_STORAGE);
    setExternalStateVariable(s, "t_esv",
                             std::span<mgis::real>(zeros.data(), 9u),
                             MaterialStateManager::EXTERNAL_STORAGE);
    setExternalStateVariable(s, "t2_esv[0]",
                             std::span<mgis::real>(zeros.data(), 9u),
                             MaterialStateManager::EXTERNAL_STORAGE);
    setExternalStateVariable(s, "t2_esv[1]",
                             std::span<mgis::real>(zeros.data(), 9u),
                             MaterialStateManager::EXTERNAL_STORAGE);
  };
  init(d.s0);
  init(d.s1);
  check(isExternalStateVariableUniform(d.s0, "Temperature"),
        "Temperature shall be uniform");
  check(isExternalStateVariableUniform(d.s0, "v_esv"),
        "v_esv shall be uniform");
  check(isExternalStateVariableUniform(d.s1, "v_esv"),
        "v_esv shall be uniform");
  check(!isExternalStateVariableUniform(d.s1, "v2_esv[0]"),
        "v2_esv[0] shall not be uniform");
  check(!isExternalStateVariableUniform(d.s1, "v2_esv[1]"),
        "v2_esv[1] shall not be uniform");
  // checking if the external storage do work as expected
  v2_esv1_values[3] = -1;
  //
  integrate(d, IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR, dt, 0, 2);
  //
  for (mgis::size_type i = 0; i != 3; ++i) {
    check(std::abs(d.s1.internal_state_variables[i] - v_esv_values.at(i)) < eps,
          "invalid internal state variable value");
    check(std::abs(d.s1.internal_state_variables[54 + i] - v_esv_values.at(i)) <
              eps,
          "invalid internal state variable value");
    check(std::abs(d.s1.internal_state_variables[i + 3] -
                   v2_esv0_values.at(i)) < eps,
          "invalid internal state variable value");
    check(std::abs(d.s1.internal_state_variables[57 + i] -
                   v2_esv0_values.at(i + 3)) < eps,
          "invalid internal state variable value");
    check(std::abs(d.s1.internal_state_variables[i + 6] -
                   v2_esv1_values.at(i)) < eps,
          "invalid internal state variable value");
    check(std::abs(d.s1.internal_state_variables[60 + i] -
                   v2_esv1_values.at(i + 3)) < eps,
          "invalid internal state variable value");
  }
}  // end of check_material_data_manager

static void check_update(const mgis::behaviour::Behaviour& b) {
  using namespace mgis::behaviour;
  using MutableFieldHolder = MaterialStateManager::MutableFieldHolder;
  constexpr const auto eps = 1.e-14;
  // s_esv given at each integration point in one state and uniform in the
  // other one
  auto uniform_values = std::vector<mgis::real>{1, 2, 3, 4, 5, 6};
  auto values = std::vector<mgis::real>(12, 0);
  MaterialDataManager d{b, 2};
  setExternalStateVariable(d.s0, "s_esv", values,
                           MaterialStateManager::EXTERNAL_STORAGE);
  setExternalStateVariable(d.s1, "s_esv", uniform_values);
  update(d);
  for (mgis::size_type i = 0; i != 12; ++i) {
    check(std::abs(values[i] - uniform_values[i % 6]) < eps,
          "invalid external state variable value after update");
  }
  setExternalStateVariable(d.s0, "s_esv", uniform_values);
  setExternalStateVariable(d.s1, "s_esv", values);
  revert(d);
  const auto& s_esv =
      std::get<MutableFieldHolder>(d.s1.external_state_variables.at("s_esv"));
  check(std::get<std::vector<mgis::real>>(s_esv.value) == uniform_values,
        "invalid external state variable value after revert");
}  // end of check_update

static void check_external_storage_hook(const mgis::behaviour::Behaviour& b) {
  using namespace mgis;
  using namespace mgis::behaviour;
  constexpr auto es = MaterialStateManager::EXTERNAL_STORAGE;
  constexpr auto ls = MaterialStateManager::LOCAL_STORAGE;
  constexpr auto update_policy = MaterialStateManager::NOUPDATE;
  const auto T = std::array<real, 2>{293.15, 293.15};
  const auto values = std::span<const real>{T};
  auto ncalls = 0;
  const auto hook = [&ncalls] { ++ncalls; };
  auto ctx = Context{};
  {
    auto m = MaterialDataManager{b, 2};
    check(setExternalStateVariable(ctx, m.s0, "Temperature", values, es,
                                   update_policy, hook) &&
              setExternalStateVariable(ctx, m.s1, "Temperature", values, es,
                                       update_policy, hook),
          "setting immutable values failed");
    update(m);
    revert(m);
    check(ncalls == 0, "immutable values released by update or revert");
    auto copy = m.s1.external_state_variables.at("Temperature");
    setExternalStateVariable(m.s1, "Temperature", 293.15);
    check(ncalls == 0,
          "immutable values released although a copy refers to them");
    copy = 293.15;
    check(ncalls == 1, "immutable values not released when replaced");
  }
  check(ncalls == 2, "immutable values not released by the destruction");
  auto m2 = MaterialDataManager{b, 2};
  check(setExternalStateVariable(ctx, m2.s1, "Temperature", values, ls,
                                 update_policy, hook),
        "setting local values failed");
  check(ncalls == 3, "hook not called at once with a local storage");
}  // end of check_external_storage_hook

int main(const int argc, const char* const* argv) {
  using namespace mgis::behaviour;
  constexpr const auto h = Hypothesis::TRIDIMENSIONAL;
  if (!check(argc == 3, "expected three arguments")) {
    return EXIT_FAILURE;
  }
  try {
    const auto b = load(argv[1], argv[2], h);
    check_behaviour(b, h);
    check_behaviour_data(b);
    check_material_data_manager(b);
    check_update(b);
    check_external_storage_hook(b);
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
