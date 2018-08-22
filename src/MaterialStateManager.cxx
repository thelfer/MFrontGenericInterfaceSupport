/*!
 * \file   MaterialStateManager.cxx
 * \brief
 * \author Thomas Helfer
 * \date   21/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/MaterialStateManager.hxx"

namespace mgis {

  namespace behaviour {

    MaterialStateManager::MaterialStateManager(MaterialStateManager&&) =
        default;
    MaterialStateManager::MaterialStateManager(const MaterialStateManager&) =
        default;

    MaterialStateManager& MaterialStateManager::operator=(
        MaterialStateManager&& src) {
      mgis::raise_if(&src.b != &this->b,
                     "MaterialStateManager::operator=: unmatched behaviour");
      if (&src != this) {
        this->gradients = std::move(src.gradients);
        this->thermodynamic_forces = std::move(src.thermodynamic_forces);
        this->material_properties = std::move(src.material_properties);
        this->internal_state_variables =
            std::move(src.internal_state_variables);
        this->external_state_variables =
            std::move(src.external_state_variables);
      }
      return *this;
    }  // end of MaterialStateManager::operator=

    MaterialStateManager& MaterialStateManager::operator=(
        const MaterialStateManager& src) {
      mgis::raise_if(&src.b != &this->b,
                     "MaterialStateManager::operator=: unmatched behaviour");
      if (&src != this) {
        this->gradients = src.gradients;
        this->thermodynamic_forces = src.thermodynamic_forces;
        this->material_properties = src.material_properties;
        this->internal_state_variables = src.internal_state_variables;
        this->external_state_variables = src.external_state_variables;
      }
      return *this;
    }  // end of MaterialStateManager::operator=

    MaterialStateManager::~MaterialStateManager() = default;

    MaterialStateManager::MaterialStateManager(const Behaviour& behaviour,
                                               const size_type s)
        : b(behaviour) {
      auto init = [this, s](std::vector<real>& values,
                            const std::vector<Variable> variables) {
        constexpr const auto zero = real{0};
        const auto as = getArraySize(variables, this->b.hypothesis);
        values.resize(s * as, zero);
      };
      init(this->gradients, this->b.gradients);
      init(this->thermodynamic_forces, this->b.thermodynamic_forces);
      //       init(this->material_properties, this->b.mps);
      init(this->internal_state_variables, this->b.isvs);
      //       init(this->external_state_variables, this->b.esvs);
    }  // end of MaterialStateManager::MaterialStateManager

  }  // end of namespace behaviour

}  // end of namespace mgis
