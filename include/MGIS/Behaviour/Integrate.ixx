/*!
 * \file   include/MGIS/Behaviour/Integrate.ixx
 * \brief    
 * \author Thomas Helfer
 * \date   01/08/2018
 * \copyright Copyright (C) 2006-2018 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence
 * or the CECILL-A licence. A copy of thoses licences are delivered
 * with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MGIS_BEHAVIOUR_INTEGRATE_IXX
#define LIB_MGIS_BEHAVIOUR_INTEGRATE_IXX

#include "MGIS/Behaviour/Behaviour.hxx"

namespace mgis{

  namespace behaviour {

    inline int integrate(BehaviourDataView& d, const Behaviour& b) {
      return b.b(&d);
    }  // end of integrate

  }  // end of namespace behaviour

} // end of namespace mgis

#endif /* LIB_MGIS_BEHAVIOUR_INTEGRATE_IXX */
