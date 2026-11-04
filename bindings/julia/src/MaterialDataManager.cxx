/*!
 * \file   bindings/julia/src/MaterialDataManager.cxx
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

#include <jlcxx/jlcxx.hpp>
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Julia/JuliaUtilities.hxx"

void declareMaterialDataManager();

void declareMaterialDataManager(jlcxx::Module& m) {
  using mgis::behaviour::MaterialDataManager;
  using mgis::behaviour::MaterialStateManager;
  m.add_type<MaterialDataManager>("MaterialDataManager")
      .method("get_s0",
              [](MaterialDataManager& d) noexcept -> MaterialStateManager& {
                return d.s0;
              })
      .method("get_s1",
              [](MaterialDataManager& d) noexcept -> MaterialStateManager& {
                return d.s1;
              })
      .method("update", [](MaterialDataManager& d) { update(d); })
      .method("revert", [](MaterialDataManager& d) { revert(d); });
  // MaterialDataManager(b, n) creates a material data manager which keeps b
  // alive
  m.add_type<mgis::julia::Wrapper<MaterialDataManager>>(
      "MaterialDataManagerWrapper",
      jlcxx::julia_base_type<MaterialDataManager>());
  m.constructor<mgis::julia::Wrapper<MaterialDataManager>,
                const mgis::julia::SharedBehaviour&, mgis::size_type>(
      jlcxx::julia_base_type<MaterialDataManager>());
}  // end of declareMaterialDataManager
