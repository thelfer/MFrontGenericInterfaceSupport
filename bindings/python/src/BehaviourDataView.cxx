/*!
 * \file   BehaviourDataView.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   08/11/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <boost/python/class.hpp>
#include "MGIS/Behaviour/BehaviourDataView.hxx"

void declareBehaviourDataView();

void declareBehaviourDataView() {
  using mgis::behaviour::BehaviourDataView;
  // exporting the BehaviourData class
  boost::python::class_<BehaviourDataView>("BehaviourDataView");
}  // end of declareBehaviourData
