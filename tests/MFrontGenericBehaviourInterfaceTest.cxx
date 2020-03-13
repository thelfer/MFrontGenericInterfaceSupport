/*!
 * \file   MFrontGenericBehaviourInterfaceTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   20/06/2018
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
    check(!isStandardFiniteStrainBehaviour(argv[1], argv[2]),
          "invalid behaviour type");
    const auto d = load(argv[1], argv[2], h);
    check(d.behaviour == "Gurson", "invalid behaviour name");
    check(d.hypothesis == h, "invalid hypothesis");
    check(d.source == "Gurson.mfront", "invalid source");
    check(d.tfel_version == TFEL_VERSION, "invalid TFEL version");
    check(d.mps.size() == 0, "invalid number of material properties");
    check(d.btype == Behaviour::STANDARDSTRAINBASEDBEHAVIOUR,
	  "invalid behaviour type");
    check(d.kinematic == Behaviour::SMALLSTRAINKINEMATIC,
	  "invalid kinematic value");
    check(d.symmetry == Behaviour::ISOTROPIC, "invalid behaviour symmetry");
    if (check(d.gradients.size() == 1u, "invalid number of gradients")) {
      check(d.gradients[0].name == "Strain", "invalid gradient name");
      check(d.gradients[0].type == Variable::STENSOR, "invalid gradient type");
    }
    if (check(d.thermodynamic_forces.size() == 1u,
	      "invalid number of thermodynamic_forces")) {
      check(d.thermodynamic_forces[0].name == "Stress", "invalid flux name");
      check(d.thermodynamic_forces[0].type == Variable::STENSOR,
	    "invalid flux type");
    }
    if (check(d.isvs.size() == 4, "invalid number of internal state variables")) {
      check(d.isvs[0].name == "ElasticStrain",
	    "invalid name for the first internal state variable");
      check(d.isvs[0].type == Variable::STENSOR,
	    "invalid type for the first internal state variable");
      check(d.isvs[1].name == "EquivalentPlasticStrain",
	    "invalid name for the second internal state variable");
      check(d.isvs[1].type == Variable::SCALAR,
	    "invalid type for the second internal state variable");
      check(d.isvs[2].name == "MatrixEquivalentPlasticStrain",
	    "invalid name for the third internal state variable");
      check(d.isvs[2].type == Variable::SCALAR,
	    "invalid type for the third internal state variable");
      check(d.isvs[3].name == "Porosity",
	    "invalid name for the third internal state variable");
      check(d.isvs[3].type == Variable::SCALAR,
	    "invalid type for the fourth internal state variable");
    }
    if (check(d.esvs.size() == 1, "invalid number of external state variables")) {
      check(d.esvs[0].name == "Temperature",
	    "invalid name for the first external state variable");
      check(d.esvs[0].type == Variable::SCALAR,
	    "invalid type for the first external state variable");
    }
  } catch(std::exception& e){
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
