/*!
 * \file   bindings/julia/src/BehaviourData.cxx
 * \brief
 * \author Thomas Helfer
 * \date   17/05/2019
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <jlcxx/jlcxx.hpp>
#include "MGIS/Behaviour/BehaviourData.hxx"
#include "MGIS/Julia/JuliaUtilities.hxx"

void declareBehaviourData();

void declareBehaviourData(jlcxx::Module& m) {
  using mgis::behaviour::BehaviourData;
  using mgis::behaviour::BehaviourDataView;
  using mgis::behaviour::State;

  m.add_type<BehaviourData>("BehaviourData")
      .constructor<const mgis::behaviour::Behaviour&>()
      .method("get_time_increment",
              [](BehaviourData& d) noexcept -> mgis::real& { return d.dt; })
      .method("set_time_increment!",
              [](BehaviourData& d, const mgis::real v) noexcept { d.dt = v; })
      .method("get_tangent_operator",
              [](BehaviourData& d) noexcept -> std::vector<mgis::real>& {
                return d.K;
              })
      .method("set_tangent_operator!",
              [](BehaviourData& d, const jlcxx::ArrayRef<mgis::real>& a) {
                mgis::julia::assign(d.K, a);
              })
      .method("get_time_increment_increase_factor",
              [](BehaviourData& d) noexcept -> mgis::real& { return d.rdt; })
      .method("set_time_increment_increase_factor!",
              [](BehaviourData& d, const mgis::real& v) noexcept -> void {
                d.rdt = v;
              })
      .method("get_s0",
              [](BehaviourData& d) noexcept -> State& { return d.s0; })
      .method("get_s1",
              [](BehaviourData& d) noexcept -> State& { return d.s1; })
      .method("get_initial_state",
              [](BehaviourData& d) noexcept -> State& { return d.s0; })
      .method("get_final_state",
              [](BehaviourData& d) noexcept -> State& { return d.s1; })
      .method("update", [](BehaviourData& d) { update(d); })
      .method("revert", [](BehaviourData& d) { revert(d); })
      .method("make_view", [](BehaviourData& d) -> BehaviourDataView {
        return make_view(d);
      });
}  // end of declareBehaviourData
