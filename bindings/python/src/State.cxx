/*!
 * \file   State.cxx
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
#include "MGIS/Behaviour/State.hxx"

void declareState(pybind11::module_&);

static pybind11::object State_getGradients(mgis::behaviour::State& s) {
  return mgis::python::wrapInNumPyArray(s.gradients);
}  // end of State_getGradients

static pybind11::object State_getThermodynamicForces(
    mgis::behaviour::State& s) {
  return mgis::python::wrapInNumPyArray(s.thermodynamic_forces);
}  // end of State_getThermodynamicForces

static pybind11::object State_getMaterialProperties(mgis::behaviour::State& s) {
  return mgis::python::wrapInNumPyArray(s.material_properties);
}  // end of State_getMaterialProperties

static pybind11::object State_getInternalStateVariables(
    mgis::behaviour::State& s) {
  return mgis::python::wrapInNumPyArray(s.internal_state_variables);
}  // end of State_getInternalStateVariables

static pybind11::object State_getExternalStateVariables(
    mgis::behaviour::State& s) {
  return mgis::python::wrapInNumPyArray(s.external_state_variables);
}  // end of State_getExternalStateVariables

static void State_setMaterialProperty(mgis::behaviour::State& s,
                                           const std::string& n,
                                           const mgis::real v) {
  mgis::behaviour::setMaterialProperty(s, n, v);
}  // end of State_setMaterialProperty

static void State_setMaterialProperty2(mgis::behaviour::State& s,
                                       const mgis::size_type o,
                                       const mgis::real v) {
  mgis::behaviour::setMaterialProperty(s, o, v);
}  // end of State_setMaterialProperty2

static void State_setExternalStateVariable(mgis::behaviour::State& s,
                                           const std::string& n,
                                           pybind11::object v) {
  if (pybind11::isinstance<pybind11::float_>(v)) {
    mgis::behaviour::setExternalStateVariable(s, n, pybind11::cast<double>(v));
  } else if (pybind11::isinstance<pybind11::int_>(v)) {
    mgis::behaviour::setExternalStateVariable(s, n, pybind11::cast<double>(v));
  } else {
    mgis::behaviour::setExternalStateVariable(
        s, n, mgis::python::mgis_convert_to_span(v));
  }
}  // end of State_setExternalStateVariable

static void State_setExternalStateVariable2(mgis::behaviour::State& s,
                                            const mgis::size_type o,
                                            pybind11::object v) {
  if (pybind11::isinstance<pybind11::float_>(v)) {
    mgis::behaviour::setExternalStateVariable(s, o, pybind11::cast<double>(v));
  } else if (pybind11::isinstance<pybind11::int_>(v)) {
    mgis::behaviour::setExternalStateVariable(s, o, pybind11::cast<double>(v));
  } else {
    mgis::behaviour::setExternalStateVariable(
        s, o, mgis::python::mgis_convert_to_span(v));
  }
}  // end of State_setExternalStateVariable2

void declareState(pybind11::module_& m) {
  using mgis::behaviour::Behaviour;
  using mgis::behaviour::State;
  pybind11::class_<State>(m, "State")
      .def_readonly("mass_density", &State::mass_density)
      .def_readonly("stored_energy", &State::stored_energy)
      .def_readonly("dissipated_energy", &State::dissipated_energy)
      .def_property_readonly("gradients", &State_getGradients)
      .def_property_readonly("thermodynamic_forces",
                             &State_getThermodynamicForces)
      .def_property_readonly("material_properties",
                             &State_getMaterialProperties)
      .def_property_readonly("internal_state_variables",
                             &State_getInternalStateVariables)
      .def_property_readonly("external_state_variables",
                             &State_getExternalStateVariables);
  //
  m.def("setMaterialProperty", State_setMaterialProperty,
        "set the value of a material property by name");
  m.def("setMaterialProperty", State_setMaterialProperty2,
        "set the value of a material property by offset");
  m.def("setExternalStateVariable", State_setExternalStateVariable,
        "set the value of an external state variable by name");
  m.def("setExternalStateVariable", State_setExternalStateVariable2,
        "set the value of an external state variable by offset");

}  // end of declareState
