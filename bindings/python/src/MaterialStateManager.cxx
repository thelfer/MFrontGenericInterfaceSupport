/*!
 * \file   MaterialStateManager.cxx
 * \brief
 * \author Thomas Helfer
 * \date   10/11/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <boost/python/def.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/class.hpp>
#include "MGIS/Raise.hxx"
#include "MGIS/Python/NumPySupport.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/MaterialStateManager.hxx"

static void MaterialStateManagerInitializer_bindGradients(
    mgis::behaviour::MaterialStateManagerInitializer& i,
    boost::python::object K) {
  i.gradients = mgis::python::mgis_convert_to_span(K);
}  // end of MaterialStateManagerInitializer_bindGradients

static void MaterialStateManagerInitializer_bindThermodynamicForces(
    mgis::behaviour::MaterialStateManagerInitializer& i,
    boost::python::object K) {
  i.thermodynamic_forces = mgis::python::mgis_convert_to_span(K);
}  // end of MaterialStateManagerInitializer_bindThermodynamicForces

static void MaterialStateManagerInitializer_bindInternalStateVariables(
    mgis::behaviour::MaterialStateManagerInitializer& i,
    boost::python::object K) {
  i.internal_state_variables = mgis::python::mgis_convert_to_span(K);
}  // end of MaterialStateManagerInitializer_bindInternalStateVariables

static void MaterialStateManagerInitializer_bindStoredEnergies(
    mgis::behaviour::MaterialStateManagerInitializer& i,
    boost::python::object K) {
  i.stored_energies = mgis::python::mgis_convert_to_span(K);
}  // end of MaterialStateManagerInitializer_bindStoredEnergies

static void MaterialStateManagerInitializer_bindDissipatedEnergies(
    mgis::behaviour::MaterialStateManagerInitializer& i,
    boost::python::object K) {
  i.dissipated_energies = mgis::python::mgis_convert_to_span(K);
}  // end of MaterialStateManagerInitializer_bindDissipatedEnergies

static boost::python::object MaterialStateManager_getGradients(
    mgis::behaviour::MaterialStateManager& s) {
  return mgis::python::wrapInNumPyArray(s.gradients, s.gradients_stride);
}  // end of MaterialStateManager_getGradients

static boost::python::object MaterialStateManager_getThermodynamicForces(
    mgis::behaviour::MaterialStateManager& s) {
  return mgis::python::wrapInNumPyArray(s.thermodynamic_forces,
                                        s.thermodynamic_forces_stride);
}  // end of MaterialStateManager_getThermodynamicForces

static boost::python::object MaterialStateManager_getInternalStateVariables(
    mgis::behaviour::MaterialStateManager& s) {
  return mgis::python::wrapInNumPyArray(s.internal_state_variables,
                                        s.internal_state_variables_stride);
}  // end of MaterialStateManager_getInternalStateVariables

static boost::python::object MaterialStateManager_getStoredEnergies(
    mgis::behaviour::MaterialStateManager& s) {
  if (!s.b.computesStoredEnergy) {
    mgis::raise(
        "MaterialStateManager_getStoredEnergies: "
        "the stored energy is not computed by the behaviour");
  }
  return mgis::python::wrapInNumPyArray(s.stored_energies);
}  // end of MaterialStateManager_getStoredEnergy

static boost::python::object MaterialStateManager_getDissipatedEnergies(
    mgis::behaviour::MaterialStateManager& s) {
  if (!s.b.computesDissipatedEnergy) {
    mgis::raise(
        "MaterialStateManager_getDissipatedEnergies: "
        "the dissipated energy is not computed by the behaviour");
  }
  return mgis::python::wrapInNumPyArray(s.dissipated_energies);
}  // end of MaterialStateManager_getDissipatedEnergy

static void MaterialStateManager_setMaterialProperty(
    mgis::behaviour::MaterialStateManager& s,
    const std::string& n,
    const mgis::real v) {
  mgis::behaviour::setMaterialProperty(s, n, v);
}  // end of MaterialStateManager_setMaterialProperty

static void MaterialStateManager_setMaterialProperty2(
    mgis::behaviour::MaterialStateManager& sm,
    const std::string& n,
    const boost::python::object& o,
    const mgis::behaviour::MaterialStateManager::StorageMode s) {
  setMaterialProperty(sm, n, mgis::python::mgis_convert_to_span(o), s);
}  // end of MaterialStateManager_setMaterialProperty

static void MaterialStateManager_setExternalStateVariable(
    mgis::behaviour::MaterialStateManager& s,
    const std::string& n,
    const mgis::real v) {
  mgis::behaviour::setExternalStateVariable(s, n, v);
}  // end of MaterialStateManager_setExternalStateVariable

static void MaterialStateManager_setExternalStateVariable2(
    mgis::behaviour::MaterialStateManager& sm,
    const std::string& n,
    const boost::python::object& o,
    const mgis::behaviour::MaterialStateManager::StorageMode s) {
  setExternalStateVariable(sm, n, mgis::python::mgis_convert_to_span(o), s);
}  // end of MaterialStateManager_setExternalStateVariable

void declareMaterialStateManager();

void declareMaterialStateManager() {
  using mgis::size_type;
  using mgis::behaviour::Behaviour;
  using mgis::behaviour::MaterialStateManager;
  using mgis::behaviour::MaterialStateManagerInitializer;
  // wrapping the MaterialStateManager::StorageMode enum
  boost::python::enum_<MaterialStateManager::StorageMode>(
      "MaterialStateManagerStorageMode")
      .value("LOCAL_STORAGE", MaterialStateManager::StorageMode::LOCAL_STORAGE)
      .value("LocalStorage", MaterialStateManager::StorageMode::LOCAL_STORAGE)
      .value("LOCALSTORAGE", MaterialStateManager::StorageMode::LOCAL_STORAGE)
      .value("EXTERNAL_STORAGE",
             MaterialStateManager::StorageMode::EXTERNAL_STORAGE)
      .value("EXTERNALSTORAGE",
             MaterialStateManager::StorageMode::EXTERNAL_STORAGE)
      .value("ExternalStorage",
             MaterialStateManager::StorageMode::EXTERNAL_STORAGE);
  // wrapping the MaterialStateManagerInitializer class
  boost::python::class_<MaterialStateManagerInitializer>(
      "MaterialStateManagerInitializer")
      .def("bindGradients", &MaterialStateManagerInitializer_bindGradients,
           "use the given array to store the gradients")
      .def("bindThermodynamicForces",
           &MaterialStateManagerInitializer_bindThermodynamicForces,
           "use the given array to store the thermodynamic forces")
      .def("bindInternalStateVariables",
           &MaterialStateManagerInitializer_bindInternalStateVariables,
           "use the given array to store the internal state variables")
      .def("bindStoredEnergies",
           &MaterialStateManagerInitializer_bindStoredEnergies,
           "use the given array to store the stored energies")
      .def("bindDissipatedEnergies",
           &MaterialStateManagerInitializer_bindDissipatedEnergies,
           "use the given array to store the dissipated energies");
  // wrapping the MaterialStateManager class
  boost::python::class_<MaterialStateManager, boost::noncopyable>(
      "MaterialStateManager",
      boost::python::init<const Behaviour&, const size_type>())
      .def(boost::python::init<const Behaviour&, const size_type,
                               const MaterialStateManagerInitializer&>())
      .def_readonly("n", &MaterialStateManager::n)
      .def_readonly("number_of_integration_points", &MaterialStateManager::n)
      .add_property("gradients", &MaterialStateManager_getGradients)
      .def_readonly("gradients_stride", &MaterialStateManager::gradients_stride)
      .add_property("thermodynamic_forces",
                    &MaterialStateManager_getThermodynamicForces)
      .def_readonly("thermodynamic_forces_stride",
                    &MaterialStateManager::thermodynamic_forces_stride)
      .add_property("stored_energies", &MaterialStateManager_getStoredEnergies)
      .add_property("dissipated_energies",
                    &MaterialStateManager_getDissipatedEnergies)
      .add_property("internal_state_variables",
                    &MaterialStateManager_getInternalStateVariables)
      .def("setMaterialProperty", &MaterialStateManager_setMaterialProperty)
      .def("setMaterialProperty", &MaterialStateManager_setMaterialProperty2)
      .def("setExternalStateVariable",
           &MaterialStateManager_setExternalStateVariable)
      .def("setExternalStateVariable",
           &MaterialStateManager_setExternalStateVariable2);
  // free functions
  boost::python::def("setMaterialProperty",
                     &MaterialStateManager_setMaterialProperty);
  boost::python::def("setMaterialProperty",
                     &MaterialStateManager_setMaterialProperty2);
  boost::python::def("setExternalStateVariable",
                     &MaterialStateManager_setExternalStateVariable);
  boost::python::def("setExternalStateVariable",
                     &MaterialStateManager_setExternalStateVariable2);

}  // end of declareMaterialStateManager
