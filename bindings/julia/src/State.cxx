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
      .method("set_stored_energy!",
              [](State& s, const mgis::real v) { s.stored_energy = v; })
      .method("get_stored_energy",
              [](State& s) -> mgis::real& { return s.stored_energy; })
      .method("set_dissipated_energy!",
              [](State& s, const mgis::real v) { s.dissipated_energy = v; })
      .method("get_dissipated_energy",
              [](State& s) -> mgis::real& { return s.dissipated_energy; })
      .method("get_gradients",
              [](State& s) -> std::vector<mgis::real>& { return s.gradients; })
      .method("set_gradients!",
              [](State& s, const jlcxx::ArrayRef<mgis::real>& a) {
                mgis::julia::assign(s.gradients, a);
              })
      .method("get_thermodynamic_forces",
              [](State& s) -> std::vector<mgis::real>& {
                return s.thermodynamic_forces;
              })
      .method("set_thermodynamic_forces!",
              [](State& s, const jlcxx::ArrayRef<mgis::real>& a) {
                mgis::julia::assign(s.thermodynamic_forces, a);
              })
      .method("get_material_properties",
              [](State& s) -> std::vector<mgis::real>& {
                return s.material_properties;
              })
      .method("set_material_properties!",
              [](State& s, const jlcxx::ArrayRef<mgis::real>& a) {
                mgis::julia::assign(s.material_properties, a);
              })
      .method("get_internal_state_variables",
              [](State& s) -> std::vector<mgis::real>& {
                return s.internal_state_variables;
              })
      .method("set_internal_state_variables!",
              [](State& s, const jlcxx::ArrayRef<mgis::real>& a) {
                mgis::julia::assign(s.internal_state_variables, a);
              })
      .method("get_external_state_variables",
              [](State& s) -> std::vector<mgis::real>& {
                return s.external_state_variables;
              })
      .method("set_external_state_variables!",
              [](State& s, const jlcxx::ArrayRef<mgis::real>& a) {
                mgis::julia::assign(s.external_state_variables, a);
              });

} // end of declareState