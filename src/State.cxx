/*!
 * \file   State.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   02/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/State.hxx"

namespace mgis {

  namespace behaviour {

    State::State() = default;
    State::State(State&&) = default;
    State::State(const State&) = default;
    State& State::operator=(State&&) = default;
    State& State::operator=(const State&) = default;

    State::State(const Behaviour& b){
      auto init = [&b](std::vector<real>& values,
                       const std::vector<Variable> variables) {
        constexpr const auto zero = real{0};
        const auto s = getArraySize(variables, b.hypothesis);
        values.resize(s, zero);
      };
      init(this->gradients, b.gradients);
      init(this->thermodynamic_forces, b.thermodynamic_forces);
      init(this->material_properties, b.mps);
      init(this->internal_state_variables, b.isvs);
      init(this->external_state_variables, b.esvs);
    }  // end of State::State

  }  // end of namespace behaviour

}  // end of namespace mgis
