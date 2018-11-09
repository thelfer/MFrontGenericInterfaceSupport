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
  using mgis::behaviour::BehaviourDataView;
  // exporting the BehaviourData class
  boost::python::class_<BehaviourData>("BehaviourData",
                                       boost::python::init<const Behaviour&>())
      .def_readwrite("dt", &BehaviourData::dt)
      .def_readwrite("rdt", &BehaviourData::rdt)
      .add_property("s0", &BehaviourData::s0)
      .add_property("s1", &BehaviourData::s1);
  // free functions
  void (*ptr_update)(BehaviourData&) = &mgis::behaviour::update;
  void (*ptr_revert)(BehaviourData&) = &mgis::behaviour::revert;
  BehaviourDataView (*ptr_make_view)(BehaviourData&) =
      &mgis::behaviour::make_view;
  boost::python::def("update", ptr_update);
  boost::python::def("revert", ptr_revert);
  boost::python::def("make_view", ptr_make_view);
}  // end of declareBehaviourData
