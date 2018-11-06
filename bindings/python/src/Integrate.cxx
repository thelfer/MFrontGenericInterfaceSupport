/*!
 * \file   Integrate.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   06/11/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <boost/python/def.hpp>
#include "MGIS/Behaviour/Integrate.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"

void declareIntegrate();

static int integrateBehaviourData(mgis::behaviour::BehaviourData& d,
                                  const mgis::behaviour::Behaviour& b) {
  auto v = mgis::behaviour::make_view(d);
  return mgis::behaviour::integrate(v, b);
}  // end of integrateBehaviourData

void declareIntegrate() {
  boost::python::def("integrate", &integrateBehaviourData);
}  // end of declareIntegrate
