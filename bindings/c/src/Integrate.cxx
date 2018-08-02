/*!
 * \file   Integrate.cxx
 * \brief    
 * \author th202608
 * \date   02/08/2018
 * \copyright Copyright (C) 2006-2018 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence
 * or the CECILL-A licence. A copy of thoses licences are delivered
 * with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include "MGIS/Behaviour/BehaviourDataView.h"
#include "MGIS/Behaviour/Behaviour.h"
#include "MGIS/Behaviour/Integrate.h"
#include "MGIS/Behaviour/Integrate.hxx"

extern "C" {

mgis_status integrate(mgis_bv_BehaviourDataView* const d,
                      const mgis_bv_Behaviour* const b) {
  const auto r = mgis::behaviour::integrate(*d, *b);
  if(r!=1){
    return mgis_report_failure("behaviour integration failed");
  }
  return mgis_report_success();
}  // end of integrate

} // end of extern "C"
