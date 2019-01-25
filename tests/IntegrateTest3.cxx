/*!
 * \file   IntegrateTest3.cxx
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
#include "MGIS/ThreadPool.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

int main(const int argc, const char* const* argv) {
  using namespace mgis;
  using namespace mgis::behaviour;
  if (argc != 2) {
    std::cerr << "IntegrateTest: invalid number of arguments\n";
    std::exit(-1);
  }
  try{
    const auto b = load(argv[1], "Norton", Hypothesis::TRIDIMENSIONAL);
    ThreadPool p{2};
    MaterialDataManager m{b, 100};
    const auto o =
      getVariableOffset(b.isvs, "EquivalentViscoplasticStrain", b.hypothesis);
    const auto de = 5.e-5;
    // initialize the external state variable
    m.s1.external_state_variables["Temperature"] = 293.15;
    // copy d.s1 in d.s0
    update(m);
    for (size_type idx = 0; idx != m.n; ++idx) {
      m.s1.gradients[idx * m.s1.gradients_stride] = de;
    }
    // integration
    auto pi = std::array<real, 21>{};  // values of the equivalent plastic strain
    // for the first integration point
    auto pe = std::array<real, 21>{};  // values of the equivalent plastic strain
    // for the last integration point
    const auto ni = size_type{o};
    const auto ne =
      size_type{(m.n - 1) * m.s0.internal_state_variables_stride + o};
    pi[0] = m.s0.internal_state_variables[ni];
    pe[0] = m.s0.internal_state_variables[ne];
    const auto dt = real(180);
    for (size_type i = 0; i != 20; ++i) {
      integrate(p, m,IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR, dt);
      update(m);
      for (size_type idx = 0; idx != m.n; ++idx) {
	m.s1.gradients[idx * m.s1.gradients_stride] += de;
      }
      pi[i + 1] = m.s1.internal_state_variables[ni];
      pe[i + 1] = m.s1.internal_state_variables[ne];
    }
    const auto p_ref = std::array<real, 21>{0,
					    1.3523277308229e-11,
					    1.0955374667213e-07,
					    5.5890770166084e-06,
					    3.2392193670428e-05,
					    6.645865307584e-05,
					    9.9676622883138e-05,
					    0.00013302758358953,
					    0.00016635821069889,
					    0.00019969195920296,
					    0.00023302522883648,
					    0.00026635857194317,
					    0.000299691903777,
					    0.0003330252373404,
					    0.00036635857063843,
					    0.00039969190397718,
					    0.00043302523730968,
					    0.00046635857064314,
					    0.00049969190397646,
					    0.00053302523730979,
					    0.00056635857064313};
    std::cout.precision(14);
    for (size_type i = 0; i != 21; ++i) {
      if (std::abs(pi[i] - p_ref[i]) > 1.e-12) {
	std::cerr << "IntegrateTest: invalid value for the equivalent "
	  "viscoplastic strain at the first integration point"
		  << "(expected '" << p_ref[i] << "', computed '" << pi[i]
		  << "')\n";
	return EXIT_FAILURE;
      }
      if (std::abs(pe[i] - p_ref[i]) > 1.e-12) {
	std::cerr << "IntegrateTest: invalid value for the equivalent "
	  "viscoplastic strain at the last integration point"
		  << "(expected '" << p_ref[i] << "', computed '" << pe[i]
		  << "')\n";
	return EXIT_FAILURE;
      }
    }
  } catch(std::exception& e){
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
