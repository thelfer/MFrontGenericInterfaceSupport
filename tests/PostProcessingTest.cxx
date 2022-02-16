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
  check(b.behaviour == "PostProcessingTest",
        "invalid behaviour name");
  check(b.hypothesis == h, "invalid hypothesis");
  check(b.source == "PostProcessingTest.mfront",
        "invalid source");
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

int main(const int argc, const char* const* argv) {
  using namespace mgis::behaviour;
  constexpr const auto h = Hypothesis::TRIDIMENSIONAL;
  if (!check(argc == 3, "expected three arguments")) {
    return EXIT_FAILURE;
  }
  try {
    const auto b = load(argv[1], argv[2], h);
    check_behaviour(b, h);
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main