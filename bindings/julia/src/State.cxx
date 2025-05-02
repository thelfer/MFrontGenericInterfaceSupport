/*!
 * \file   bindings/julia/src/State.cxx
 * \brief
 * \author Thomas Helfer
 * \date   16/05/2019
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <jlcxx/jlcxx.hpp>
#include "MGIS/Behaviour/State.hxx"
#include "MGIS/Julia/JuliaUtilities.hxx"

void declareState();

void declareState(jlcxx::Module& m) {
  using mgis::behaviour::Behaviour;
  using mgis::behaviour::State;
  m.add_type<State>("State")
      .constructor<const Behaviour&>()
      .method(
          "set_stored_energy!",
          [](State& s, const mgis::real v) noexcept { s.stored_energy = v; })
      .method("get_stored_energy",
              [](State& s) noexcept -> mgis::real& { return s.stored_energy; })
      .method("set_dissipated_energy!",
              [](State& s, const mgis::real v) noexcept {
                s.dissipated_energy = v;
              })
      .method(
          "get_dissipated_energy",
          [](State& s) noexcept -> mgis::real& { return s.dissipated_energy; })
      .method("get_gradients",
              [](State& s) noexcept -> std::vector<mgis::real>& {
                return s.gradients;
              })
      .method("set_gradients!",
              [](State& s, const jlcxx::ArrayRef<mgis::real>& a) {
                mgis::julia::assign(s.gradients, a);
              })
      .method("get_thermodynamic_forces",
              [](State& s) noexcept -> std::vector<mgis::real>& {
                return s.thermodynamic_forces;
              })
      .method("set_thermodynamic_forces!",
              [](State& s, const jlcxx::ArrayRef<mgis::real>& a) {
                mgis::julia::assign(s.thermodynamic_forces, a);
              })
      .method("get_material_properties",
              [](State& s) noexcept -> std::vector<mgis::real>& {
                return s.material_properties;
              })
      .method("set_material_properties!",
              [](State& s, const jlcxx::ArrayRef<mgis::real>& a) {
                mgis::julia::assign(s.material_properties, a);
              })
      .method("get_internal_state_variables",
              [](State& s) noexcept -> std::vector<mgis::real>& {
                return s.internal_state_variables;
              })
      .method("set_internal_state_variables!",
              [](State& s, const jlcxx::ArrayRef<mgis::real>& a) {
                mgis::julia::assign(s.internal_state_variables, a);
              })
      .method("get_external_state_variables",
              [](State& s) noexcept -> std::vector<mgis::real>& {
                return s.external_state_variables;
              })
      .method("set_external_state_variables!",
              [](State& s, const jlcxx::ArrayRef<mgis::real>& a) {
                mgis::julia::assign(s.external_state_variables, a);
              })
      .method("set_gradient!",
              [](State& s, const std::string& n, const mgis::real v) {
                mgis::behaviour::setGradient(s, n, v);
              })
      .method("set_thermodynamic_force!",
              [](State& s, const std::string& n, const mgis::real v) {
                mgis::behaviour::setThermodynamicForce(s, n, v);
              })
      .method("set_material_property!",
              [](State& s, const std::string& n, const mgis::real v) {
                mgis::behaviour::setMaterialProperty(s, n, v);
              })
      .method("set_internal_state_variable!",
              [](State& s, const std::string& n, const mgis::real v) {
                mgis::behaviour::setInternalStateVariable(s, n, v);
              })
      .method("set_external_state_variable!",
              [](State& s, const std::string& n, const mgis::real v) {
                mgis::behaviour::setExternalStateVariable(s, n, v);
              });

}  // end of declareState
