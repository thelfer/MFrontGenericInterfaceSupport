/*!
 * \file   tests/IntegrateTest4.cxx
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
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/ThreadPool.hxx"
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
    ThreadPool p{2};
    MaterialDataManager m{model, 100};
    const auto o = getVariableOffset(model.isvs, "x", model.hypothesis);
    // initialize the internal state variable
    auto pos = o;
    auto& isvs = m.s1.internal_state_variables;
    for (size_type i = 0; i != 100; ++i) {
      isvs[pos] = 1;
      pos += m.s1.internal_state_variables_stride;
    }
    // initialize the external state variable
    m.s1.external_state_variables["Temperature"] = 293.15;
    // copy d.s1 in d.s0
    update(m);
    // values of x for the first integration point
    auto xvi = std::array<real, 11>{};
    // values of x for the last integration point
    auto xve = std::array<real, 11>{};
    const auto ni = size_type{o};
    const auto ne =
        size_type{(m.n - 1) * m.s0.internal_state_variables_stride + o};
    xvi[0] = m.s0.internal_state_variables[ni];
    xve[0] = m.s0.internal_state_variables[ne];
    const auto dt = real(0.1);
    for (size_type i = 0; i != 10; ++i) {
      integrate(p, m, IntegrationType::INTEGRATION_NO_TANGENT_OPERATOR, dt);
      update(m);
      xvi[i + 1] = m.s1.internal_state_variables[ni];
      xve[i + 1] = m.s1.internal_state_variables[ne];
    }
    std::cerr.precision(14);
    auto t = mgis::real{};
    for (size_type i = 0; i != 10; ++i) {
      constexpr auto eps = mgis::real{1e-10};
      const auto x_ref = exp(-A * t);
      if (std::abs(xvi[i] - x_ref) > eps) {
        std::cerr << "IntegrateTest: invalid value for x "
                  << " at the first integration point"
                  << "(expected '" << x_ref << "', computed '" << xvi[i]
                  << "')\n";
        return EXIT_FAILURE;
      }
      if (std::abs(xve[i] - x_ref) > eps) {
        std::cerr << "IntegrateTest: invalid value for x "
                  << "at the last integration point"
                  << "(expected '" << x_ref << "', computed '" << xve[i]
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
