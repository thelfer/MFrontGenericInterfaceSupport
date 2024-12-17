/*!
 * \file   tests/IntegrateTest5.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/11/2021
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
#include "MGIS/Behaviour/BehaviourData.hxx"
#include "MGIS/Model/Model.hxx"
#include "MGIS/Behaviour/Integrate.hxx"

int main(const int argc, const char* const* argv) {
  using namespace mgis;
  using namespace mgis::behaviour;
  if (argc != 2) {
    std::cerr << "IntegrateTest: invalid number of arguments\n";
    std::exit(-1);
  }
  try {
    const auto model =
        mgis::model::load(argv[1], "ode_rk54", Hypothesis::TRIDIMENSIONAL);
    const auto A = getParameterDefaultValue<mgis::real>(model, "A");
    auto d = BehaviourData{model};
    const auto o = getVariableOffset(model.isvs, "x", model.hypothesis);
    // initialize the internal state variable
    d.s1.internal_state_variables[o] = 1;
    // initialize the external state variable
    setExternalStateVariable(d.s1, "Temperature", 293.15);
    // copy d.s1 in d.s0
    update(d);
    // values of 'x' at each time step
    auto xvalues = std::array<real, 11>{};
    xvalues[0] = d.s0.internal_state_variables[o];
    const auto dt = real(0.1);
    d.K[0] = 0;
    for (size_type i = 0; i != 10; ++i) {
      d.dt = dt;
      auto v = make_view(d);
      if (integrate(v, model) != 1) {
        return EXIT_FAILURE;
      }
      update(d);
      xvalues[i + 1] = d.s1.internal_state_variables[o];
    }
    std::cerr.precision(14);
    auto t = mgis::real{};
    for (size_type i = 0; i != 10; ++i) {
      constexpr auto eps = mgis::real{1e-10};
      const auto x_ref = exp(-A * t);
      if (std::abs(xvalues[i] - x_ref) > eps) {
        std::cerr << "IntegrateTest: invalid value for x "
                  << "at the first integration point"
                  << "(expected '" << x_ref << "', computed '" << xvalues[i]
                  << "')\n";
        return EXIT_FAILURE;
      }
      t += dt;
    }
  } catch (std::exception& e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
