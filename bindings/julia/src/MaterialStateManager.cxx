/*!
 * \file   bindings/julia/src/MaterialStateManager.cxx
 * \brief
 * \author Thomas Helfer
 * \date   19/05/2019
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <span>
#include <jlcxx/jlcxx.hpp>
#include "MGIS/Behaviour/MaterialStateManager.hxx"
#include "MGIS/Julia/JuliaUtilities.hxx"

void declareMaterialStateManager();

//! \brief Julia array sharing the values of the given span
static jlcxx::ArrayRef<mgis::real> make_array(std::span<mgis::real> v) {
  return jlcxx::ArrayRef<mgis::real>(v.data(), v.size());
}  // end of make_array

void declareMaterialStateManager(jlcxx::Module& m) {
  using mgis::behaviour::MaterialStateManager;
  m.add_type<MaterialStateManager>("MaterialStateManager")
      .method("get_gradients",
              [](MaterialStateManager& s) { return make_array(s.gradients); })
      .method("get_thermodynamic_forces",
              [](MaterialStateManager& s) {
                return make_array(s.thermodynamic_forces);
              })
      .method("get_internal_state_variables",
              [](MaterialStateManager& s) {
                return make_array(s.internal_state_variables);
              })
      .method("set_material_property!",
              [](MaterialStateManager& s, const std::string& n,
                 const mgis::real v) {
                mgis::behaviour::setMaterialProperty(s, n, v);
              })
      .method("set_external_state_variable!",
              [](MaterialStateManager& s, const std::string& n,
                 const mgis::real v) {
                mgis::behaviour::setExternalStateVariable(s, n, v);
              });
}  // end of declareMaterialStateManager
