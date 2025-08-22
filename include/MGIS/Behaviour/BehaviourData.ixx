/*!
 * \file   BehaviourData.ixx
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

#ifndef LIB_MGIS_BEHAVIOUR_BEHAVIOURDATA_IXX
#define LIB_MGIS_BEHAVIOUR_BEHAVIOURDATA_IXX

namespace mgis::behaviour {

  inline BehaviourDataView make_view(BehaviourData& d) {
    auto get_ptr = [](std::vector<real>& v) -> real* {
      if (v.empty()) {
        return nullptr;
      }
      return &v[0];
    };  // end of get_ptr
    BehaviourDataView v;
    v.error_message = d.error_message;
    v.dt = d.dt;
    v.rdt = &(d.rdt);
    v.speed_of_sound = &(d.speed_of_sound);
    v.K = get_ptr(d.K);
    v.s0 = make_view(static_cast<const State&>(d.s0));
    v.s1 = make_view(d.s1);
    return v;
  }  // end of make_view

}  // end of namespace mgis::behaviour


#endif /* LIB_MGIS_BEHAVIOUR_BEHAVIOURDATA_IXX */
