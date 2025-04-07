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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "MGIS/Python/NumPySupport.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"

void declareBehaviourData(pybind11::module_& m);

static pybind11::object BehaviourData_getK(
    mgis::behaviour::BehaviourData& d) {
  if (d.s0.b.to_blocks.size() == 1u) {
    const auto s =
        getVariableSize(d.s0.b.to_blocks.front().first, d.s0.b.hypothesis);
    return mgis::python::wrapInNumPyArray(d.K, s);
  }
  return mgis::python::wrapInNumPyArray(d.K);
}  // end of MaterialStateManager_getK

void declareBehaviourData(pybind11::module_& m) {
  using mgis::behaviour::Behaviour;
  using mgis::behaviour::BehaviourData;
  using mgis::behaviour::BehaviourDataView;
  // pointers to free functions to disambiguate the function resolution
  void (*ptr_update)(BehaviourData&) = &mgis::behaviour::update;
  void (*ptr_revert)(BehaviourData&) = &mgis::behaviour::revert;
  // exporting the BehaviourData class
  pybind11::class_<BehaviourData>(m, "BehaviourData")
      .def(pybind11::init<const Behaviour&>())
      .def_readwrite("dt", &BehaviourData::dt)
      .def_readwrite("rdt", &BehaviourData::rdt)
      .def_readonly("s0", &BehaviourData::s0)
      .def_readonly("s1", &BehaviourData::s1)
      .def_property_readonly("K", &BehaviourData_getK)
      .def("update", ptr_update)
      .def("revert", ptr_revert);
  // free functions
  BehaviourDataView (*ptr_make_view)(BehaviourData&) =
      &mgis::behaviour::make_view;
  m.def("update", ptr_update);
  m.def("revert", ptr_revert);
  m.def("make_view", ptr_make_view);
}  // end of declareBehaviourData
