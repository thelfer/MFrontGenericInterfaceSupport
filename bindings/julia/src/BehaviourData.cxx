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

#include <memory>
#include <jlcxx/jlcxx.hpp>
#include "MGIS/Behaviour/BehaviourData.hxx"
#include "MGIS/Julia/JuliaUtilities.hxx"

namespace mgis::julia {

  /*!
   * \brief behaviour data created from Julia: they own a copy of their
   * behaviour, which is thus not destroyed if the Julia object of the
   * behaviour is collected.
   */
  struct OwningBehaviourData : mgis::behaviour::BehaviourData {
    OwningBehaviourData(const mgis::behaviour::Behaviour& b)
        : OwningBehaviourData(
              std::make_shared<const mgis::behaviour::Behaviour>(b)) {}
    OwningBehaviourData(std::shared_ptr<const mgis::behaviour::Behaviour> p)
        : mgis::behaviour::BehaviourData(*p), behaviour(std::move(p)) {}
    //! \brief behaviour, shared with the copies of this object
    std::shared_ptr<const mgis::behaviour::Behaviour> behaviour;
  };

}  // end of namespace mgis::julia

namespace jlcxx {

  //! \brief
  template <>
  struct SuperType<mgis::julia::OwningBehaviourData> {
    using type = mgis::behaviour::BehaviourData;
  };

}  // end of namespace jlcxx

void declareBehaviourData();

void declareBehaviourData(jlcxx::Module& m) {
  using mgis::behaviour::BehaviourData;
  using mgis::behaviour::BehaviourDataView;
  using mgis::behaviour::State;

  m.add_type<BehaviourData>("BehaviourData")
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
  // BehaviourData(b) creates behaviour data owning a copy of b
  m.add_type<mgis::julia::OwningBehaviourData>(
      "OwningBehaviourData", jlcxx::julia_base_type<BehaviourData>());
  m.constructor<mgis::julia::OwningBehaviourData,
                const mgis::behaviour::Behaviour&>(
      jlcxx::julia_base_type<BehaviourData>());
}  // end of declareBehaviourData
