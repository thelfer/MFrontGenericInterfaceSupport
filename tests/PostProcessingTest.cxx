/*!
 * \file   tests/PostProcessingTest.cxx
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
#include <array>
#include <cstdlib>
#include <iostream>
#include <string_view>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

bool success = true;

bool check(const bool b, const std::string_view msg) {
  if (!b) {
    success = false;
    std::cerr << msg << '\n';
  }
  return b;
}

void check_behaviour(const mgis::behaviour::Behaviour& b,
                     const mgis::behaviour::Hypothesis h) {
  using namespace mgis::behaviour;
  check(b.behaviour == "PostProcessingTest", "invalid behaviour name");
  check(b.hypothesis == h, "invalid hypothesis");
  check(b.source == "PostProcessingTest.mfront", "invalid source");
  check(b.tfel_version == TFEL_VERSION, "invalid TFEL version");
  check(b.postprocessings.size() == 1u, "invalid number of post-processings");
  check(b.postprocessings.find("PrincipalStrain") != b.postprocessings.end(),
        "invalid post-processing name");
  const auto& p = b.postprocessings.at("PrincipalStrain");
  check(p.outputs.size() == 1u, "invalid post-processing output size");
  const auto& o = p.outputs.at(0);
  check(o.name == "PrincipalStrain", "invalid post-processing output");
  check(o.type == Variable::VECTOR_3D, "invalid post-processing output type");
  check(o.type_identifier == 2 + (3 << 3),
        "invalid post-processing output type identifier");
  check(getVariableSize(o, h) == 3, "invalid post-processing output size");
}  // end of check_behaviour

void call_postprocessing(const mgis::behaviour::Behaviour& b) {
  using namespace mgis::behaviour;
  constexpr auto e =
      std::array<mgis::real, 6u>{1.3e-2, 1.2e-2, 1.4e-2, 0., 0., 0.};
  constexpr auto e2 =
      std::array<mgis::real, 6u>{1.2e-2, 1.3e-2, 1.4e-2, 0., 0., 0.};
  constexpr auto eps = 10 * std::numeric_limits<mgis::real>::epsilon();
  auto d = BehaviourData{b};
  // initialize the states
  setExternalStateVariable(d.s0, "Temperature", 293.15);
  setExternalStateVariable(d.s1, "Temperature", 293.15);
  //
  for (mgis::size_type i = 0; i != 6; ++i) {
    d.s1.gradients[i] = e[i];
  }
  //
  auto outputs = allocatePostProcessingVariables(b, "PrincipalStrain");
  if (!check(outputs.size() == 3u, "invalid outputs initialisation")) {
    return;
  }
  auto v = make_view(d);
  executePostProcessing(outputs, v, b, "PrincipalStrain");
  for (mgis::size_type i = 0; i != 3; ++i) {
    check(std::abs(outputs[i] - e2[i]) < eps, "invalid output value");
  }
}  // end of call_postprocessing

void call_postprocessing2(const mgis::behaviour::Behaviour& b) {
  using namespace mgis::behaviour;
  constexpr auto e =
      std::array<mgis::real, 6u>{1.3e-2, 1.2e-2, 1.4e-2, 0., 0., 0.};
  constexpr auto e2 =
      std::array<mgis::real, 6u>{1.2e-2, 1.3e-2, 1.4e-2, 0., 0., 0.};
  constexpr auto eps = 10 * std::numeric_limits<mgis::real>::epsilon();
  auto m = MaterialDataManager{b, 2u};
  // initialize the states
  setMaterialProperty(m.s1, "YoungModulus", 150e9);
  setMaterialProperty(m.s1, "PoissonRatio", 0.3);
  setExternalStateVariable(m.s1, "Temperature", 293.15);
  update(m);
  //
  for (mgis::size_type i = 0; i != 6; ++i) {
    m.s1.gradients[i] = e[i];
    m.s1.gradients[6 + i] = e[i];
  }
  //
  auto outputs = allocatePostProcessingVariables(m, "PrincipalStrain");
  if (!check(outputs.size() == 6u, "invalid outputs initialisation")) {
    return;
  }
  executePostProcessing(outputs, m, "PrincipalStrain");
  for (mgis::size_type i = 0; i != 3; ++i) {
    check(std::abs(outputs[i] - e2[i]) < eps, "invalid output value");
    check(std::abs(outputs[3 + i] - e2[i]) < eps, "invalid output value");
  }
}  // end of call_postprocessing2

int main(const int argc, const char* const* argv) {
  using namespace mgis::behaviour;
  constexpr const auto h = Hypothesis::TRIDIMENSIONAL;
  if (!check(argc == 3, "expected three arguments")) {
    return EXIT_FAILURE;
  }
  try {
    const auto b = load(argv[1], argv[2], h);
    check_behaviour(b, h);
    call_postprocessing(b);
    call_postprocessing2(b);
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main