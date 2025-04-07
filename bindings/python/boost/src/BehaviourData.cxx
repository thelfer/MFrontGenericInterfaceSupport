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
#include "MGIS/Python/NumPySupport.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"

void declareBehaviourData();

static boost::python::object BehaviourData_getK(
    mgis::behaviour::BehaviourData& d) {
  if (d.s0.b.to_blocks.size() == 1u) {
    const auto s =
        getVariableSize(d.s0.b.to_blocks.front().first, d.s0.b.hypothesis);
    return mgis::python::wrapInNumPyArray(d.K, s);
  }
  return mgis::python::wrapInNumPyArray(d.K);
}  // end of MaterialStateManager_getK

void declareBehaviourData() {
  using mgis::behaviour::Behaviour;
  using mgis::behaviour::BehaviourData;
  using mgis::behaviour::BehaviourDataView;
  // pointers to free functions to disambiguate the function resolution
  void (*ptr_update)(BehaviourData&) = &mgis::behaviour::update;
  void (*ptr_revert)(BehaviourData&) = &mgis::behaviour::revert;
  // exporting the BehaviourData class
  boost::python::class_<BehaviourData>("BehaviourData",
                                       boost::python::init<const Behaviour&>())
      .def_readwrite("dt", &BehaviourData::dt)
      .def_readwrite("rdt", &BehaviourData::rdt)
      .add_property("s0", &BehaviourData::s0)
      .add_property("s1", &BehaviourData::s1)
      .add_property("K", &BehaviourData_getK)
      .def("update", ptr_update)
      .def("revert", ptr_revert);
  // free functions
  BehaviourDataView (*ptr_make_view)(BehaviourData&) =
      &mgis::behaviour::make_view;
  boost::python::def("update", ptr_update);
  boost::python::def("revert", ptr_revert);
  boost::python::def("make_view", ptr_make_view);
}  // end of declareBehaviourData
