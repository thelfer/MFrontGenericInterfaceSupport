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
  m.add_type<BehaviourData>("BehaviourData")
      .constructor<const mgis::behaviour::Behaviour&>()
      .method("get_time_increment",
              [](BehaviourData& d) -> mgis::real& { return d.dt; })
      .method("set_time_increment!",
              [](BehaviourData& d, const mgis::real v) { d.dt = v; })
      .method("get_tangent_operator",
              [](BehaviourData& d) -> std::vector<mgis::real>& { return d.K; });
  //       .method("set_tangent_operator!",
  //               [](BehaviourData& d, const std::vector<mgis::real>& K) {
  //                 if (d.K.size() != K.size()) {
  //                   mgis::raise<std::range_error>(
  //                       "set_tangent_operator!: "
  //                       "unmatched size");
  //                 }
  //                 std::copy(K.begin(), K.end(), d.K.begin());
  //               })
  //       .method("set_tangent_operator!",
  //               [](BehaviourData& d, const jlcxx::ArrayRef<mgis::real>& K) {
  //                 if (d.K.size() != K.size()) {
  //                   mgis::raise<std::range_error>(
  //                       "set_tangent_operator!: "
  //                       "unmatched size");
  //                 }
  //                 std::copy(K.begin(), K.end(), d.K.begin());
  //               });
}  // end of declareBehaviourData