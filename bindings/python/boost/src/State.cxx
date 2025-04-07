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

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include "MGIS/Python/NumPySupport.hxx"
#include "MGIS/Behaviour/State.hxx"

void declareState();

static boost::python::object State_getGradients(mgis::behaviour::State& s) {
  return mgis::python::wrapInNumPyArray(s.gradients);
}  // end of State_getGradients

static boost::python::object State_getThermodynamicForces(
    mgis::behaviour::State& s) {
  return mgis::python::wrapInNumPyArray(s.thermodynamic_forces);
}  // end of State_getThermodynamicForces

static boost::python::object State_getMaterialProperties(
    mgis::behaviour::State& s) {
  return mgis::python::wrapInNumPyArray(s.material_properties);
}  // end of State_getMaterialProperties

static boost::python::object State_getInternalStateVariables(
    mgis::behaviour::State& s) {
  return mgis::python::wrapInNumPyArray(s.internal_state_variables);
}  // end of State_getInternalStateVariables

static boost::python::object State_getExternalStateVariables(
    mgis::behaviour::State& s) {
  return mgis::python::wrapInNumPyArray(s.external_state_variables);
}  // end of State_getExternalStateVariables

static void State_setExternalStateVariable(mgis::behaviour::State& s,
                                           const std::string& n,
                                           boost::python::object v) {
  auto e = boost::python::extract<double>(v);
  if (e.check()) {
    mgis::behaviour::setExternalStateVariable(s, n, e());
  } else {
    mgis::behaviour::setExternalStateVariable(
        s, n, mgis::python::mgis_convert_to_span(v));
  }
}  // end of State_setExternalStateVariable

static void State_setExternalStateVariable2(mgis::behaviour::State& s,
                                            const mgis::size_type o,
                                            boost::python::object v) {
  auto e = boost::python::extract<double>(v);
  if (e.check()) {
    mgis::behaviour::setExternalStateVariable(s, o, e());
  } else {
    mgis::behaviour::setExternalStateVariable(
        s, o, mgis::python::mgis_convert_to_span(v));
  }
}  // end of State_setExternalStateVariable2

void declareState() {
  using mgis::behaviour::Behaviour;
  using mgis::behaviour::State;
  boost::python::class_<State>("State", boost::python::no_init)
      .add_property("mass_density", &State::mass_density)
      .add_property("stored_energy", &State::stored_energy)
      .add_property("dissipated_energy", &State::dissipated_energy)
      .add_property("gradients", &State_getGradients)
      .add_property("thermodynamic_forces", &State_getThermodynamicForces)
      .add_property("material_properties", &State_getMaterialProperties)
      .add_property("internal_state_variables",
                    &State_getInternalStateVariables)
      .add_property("external_state_variables",
                    &State_getExternalStateVariables);
  //
  boost::python::def("setExternalStateVariable", State_setExternalStateVariable,
                     "set the value of an external state variable by name");
  boost::python::def("setExternalStateVariable",
                     State_setExternalStateVariable2,
                     "set the value of an external state variable by offset");

}  // end of declareState