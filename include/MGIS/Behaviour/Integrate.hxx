/*!
 * \file   Integrate.hxx
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

#ifndef LIB_MGIS_BEHAVIOUR_INTEGRATE_HXX
#define LIB_MGIS_BEHAVIOUR_INTEGRATE_HXX

#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourDataView.hxx"

namespace mgis{
  
  namespace behaviour{

    /*!
     * \brief integrate the behaviour
     * \param[in,out] d: behaviour data
     * \param[in,out] b: behaviour
     */
    int integrate(BehaviourDataView&, const Behaviour&);

  } // end of namespace behaviour

} // end of namespace mgis

#include "MGIS/Behaviour/Integrate.ixx"

#endif /* LIB_MGIS_BEHAVIOUR_INTEGRATE_HXX */
