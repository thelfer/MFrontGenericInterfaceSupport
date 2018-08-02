/*!
 * \file   BehaviourDataView.cxx
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

#include "MGIS/Behaviour/BehaviourData.hxx"
#include "MGIS/Behaviour/BehaviourDataView.hxx"

namespace mgis {

  namespace behaviour {

    BehaviourDataView make_view(BehaviourData& d) {
      auto get_ptr = [](const std::vector<real>& v) -> real* {
        if (v.empty()) {
          return nullptr;
        }
        return &v[0];
      };  // end of get_ptr
      BehaviourDataView v;
      v.dt = d.dt;
      v.rdt = d.rdt;
      v.K = get_ptr(d.K);
      v.s0 = make_view(d.s0);
      v.s1 = make_view(d.s1);
      return v;
    }  // end of make_view

  }  // end of namepace behaviour

}  // end of namespace mgis