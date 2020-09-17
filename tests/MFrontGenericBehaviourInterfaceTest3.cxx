/*!
 * \file   MFrontGenericBehaviourInterfaceTest3.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/05/2019
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <iostream>
#include <cstdlib>
#include "MGIS/Behaviour/Behaviour.hxx"

int main(const int argc, const char* const* argv) {
  using namespace mgis::behaviour;
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
  try{
    auto o = FiniteStrainBehaviourOptions{};
    o.stress_measure = FiniteStrainBehaviourOptions::PK1;
    o.tangent_operator = FiniteStrainBehaviourOptions::DPK1_DF;
    const auto d = load(o, argv[1], argv[2], h);
    check(d.symmetry == Behaviour::ORTHOTROPIC, "invalid behaviour symmetry");
    if (check(d.gradients.size() == 1u, "invalid number of gradients")) {
      check(d.gradients[0].name == "DeformationGradient", "invalid gradient name");
      check(d.gradients[0].type == Variable::TENSOR, "invalid gradient type");
    }
    if (check(d.thermodynamic_forces.size() == 1u,
	      "invalid number of thermodynamic_forces")) {
      check(d.thermodynamic_forces[0].name == "FirstPiolaKirchhoffStress",
	    "invalid flux name");
      check(d.thermodynamic_forces[0].type == Variable::TENSOR,
	    "invalid flux type");
    }
  } catch(std::exception& e){
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
