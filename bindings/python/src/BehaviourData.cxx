/*!
 * \file   bindings/python/src/BehaviourData.cxx
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
#include <boost/python/class.hpp>
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"

void declareBehaviourData();

void declareBehaviourData() {
  using mgis::behaviour::Behaviour;
  using mgis::behaviour::BehaviourData;
  // exporting the BehaviourData class
  boost::python::class_<BehaviourData>("BehaviourData",
                                       boost::python::init<const Behaviour&>())
      .add_property("dt", &BehaviourData::dt)
      .add_property("rdt", &BehaviourData::rdt)
      .add_property("s0", &BehaviourData::s0)
      .add_property("s1", &BehaviourData::s1);
  // free functions
  void (*ptr_update)(BehaviourData&) = &mgis::behaviour::update;
  void (*ptr_revert)(BehaviourData&) = &mgis::behaviour::revert;
  boost::python::def("update", ptr_update);
  boost::python::def("revert", ptr_revert);
}  // end of declareBehaviourData
