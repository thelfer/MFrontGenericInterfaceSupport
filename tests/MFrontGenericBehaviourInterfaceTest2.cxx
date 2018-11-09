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
    const auto d = load(argv[1], argv[2], h);
    check(d.behaviour == "FiniteStrainSingleCrystal", "invalid behaviour name");
    check(d.hypothesis == h, "invalid hypothesis");
    check(d.source == "FiniteStrainSingleCrystal.mfront", "invalid source");
    check(d.tfel_version == TFEL_VERSION, "invalid TFEL version");
    check(d.btype == Behaviour::STANDARDFINITESTRAINBEHAVIOUR,
	  "invalid behaviour type");
    check(d.kinematic == Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY,
	  "invalid kinematic value");
    check(d.symmetry == Behaviour::ORTHOTROPIC, "invalid behaviour symmetry");
    if (check(d.gradients.size() == 1u, "invalid number of gradients")) {
      check(d.gradients[0].name == "DeformationGradient", "invalid gradient name");
      check(d.gradients[0].type == Variable::TENSOR, "invalid gradient type");
    }
    if (check(d.thermodynamic_forces.size() == 1u,
	      "invalid number of thermodynamic_forces")) {
      check(d.thermodynamic_forces[0].name == "Stress", "invalid flux name");
      check(d.thermodynamic_forces[0].type == Variable::STENSOR,
	    "invalid flux type");
    }
    if (check(d.mps.size() == 16, "invalid number of material properties")) {
      auto check_mp = [&check](const Variable& mp, const std::string& n) {
	check(mp.name == n,
	      "invalid material property name, expected '" + n + "'");
      };
      for (const auto& mp : d.mps) {
	check(mp.type == Variable::SCALAR,
	      "invalid material property type '" + mp.name + "'");
      }
      check_mp(d.mps[0],"YoungModulus1");
      check_mp(d.mps[1],"YoungModulus2");
      check_mp(d.mps[2],"YoungModulus3");
      check_mp(d.mps[3],"PoissonRatio12");
      check_mp(d.mps[4],"PoissonRatio23");
      check_mp(d.mps[5],"PoissonRatio13");
      check_mp(d.mps[6],"ShearModulus12");
      check_mp(d.mps[7],"ShearModulus23");
      check_mp(d.mps[8], "ShearModulus13");
      check_mp(d.mps[9], "m");
      check_mp(d.mps[10], "K");
      check_mp(d.mps[11], "C");
      check_mp(d.mps[12], "R0");
      check_mp(d.mps[13], "Q");
      check_mp(d.mps[14], "b");
      check_mp(d.mps[15], "d1");
    }
    if (check(d.isvs.size() == 37, "invalid number of internal state variables")) {
      auto check_isv = [&check](const Variable& isv, const std::string& n,
				const unsigned short i) {
	const auto vn = n + "[" + std::to_string(i) + "]";
	check(isv.name == vn,
	      "invalid internal state variable name, expected '" + vn + "'");
	check(isv.type == Variable::SCALAR,
	      "invalid type for internal state variable '" + vn + "'");
      };
      for (unsigned short i = 0; i != 12; ++i) {
	check_isv(d.isvs[i], "g", i);
	check_isv(d.isvs[13 + i], "p", i);
	check_isv(d.isvs[25 + i], "a", i);
      }
      check(d.isvs[12].name == "Fe",
	    "invalid name for internal state variable 'Fe'");
      check(d.isvs[12].type == Variable::TENSOR,
	    "invalid type for internal state variable 'Fe'");
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
