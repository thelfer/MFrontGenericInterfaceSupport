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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "MGIS/Raise.hxx"
#include "MGIS/Python/NumPySupport.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/MaterialStateManager.hxx"

static void MaterialStateManagerInitializer_bindGradients(
    pybind11::object self, pybind11::object K) {
  auto& i = self.cast<mgis::behaviour::MaterialStateManagerInitializer&>();
  i.gradients = mgis::python::mgis_convert_to_bound_span(K, self);
}  // end of MaterialStateManagerInitializer_bindGradients

static void MaterialStateManagerInitializer_bindThermodynamicForces(
    pybind11::object self, pybind11::object K) {
  auto& i = self.cast<mgis::behaviour::MaterialStateManagerInitializer&>();
  i.thermodynamic_forces = mgis::python::mgis_convert_to_bound_span(K, self);
}  // end of MaterialStateManagerInitializer_bindThermodynamicForces

static void MaterialStateManagerInitializer_bindInternalStateVariables(
    pybind11::object self, pybind11::object K) {
  auto& i = self.cast<mgis::behaviour::MaterialStateManagerInitializer&>();
  i.internal_state_variables =
      mgis::python::mgis_convert_to_bound_span(K, self);
}  // end of MaterialStateManagerInitializer_bindInternalStateVariables

static void MaterialStateManagerInitializer_bindStoredEnergies(
    pybind11::object self, pybind11::object K) {
  auto& i = self.cast<mgis::behaviour::MaterialStateManagerInitializer&>();
  i.stored_energies = mgis::python::mgis_convert_to_bound_span(K, self);
}  // end of MaterialStateManagerInitializer_bindStoredEnergies

static void MaterialStateManagerInitializer_bindDissipatedEnergies(
    pybind11::object self, pybind11::object K) {
  auto& i = self.cast<mgis::behaviour::MaterialStateManagerInitializer&>();
  i.dissipated_energies = mgis::python::mgis_convert_to_bound_span(K, self);
}  // end of MaterialStateManagerInitializer_bindDissipatedEnergies

/*!
 * \brief wrap values stored per integration point
 * \param[in] values: values
 * \param[in] n: number of integration points
 * \param[in] stride: number of values per integration point
 * \param[in] base: object owning the values, kept alive by the ndarray
 */
static pybind11::object MaterialStateManager_wrapPerIntegrationPointValues(
    std::span<mgis::real>& values,
    const mgis::size_type n,
    const mgis::size_type stride,
    pybind11::handle base) {
  if (stride == 0) {
    // one empty row per integration point
    return pybind11::array_t<double>(
        std::vector<pybind11::ssize_t>{static_cast<pybind11::ssize_t>(n), 0});
  }
  return mgis::python::wrapInNumPyArray(values, stride, base);
}  // end of MaterialStateManager_wrapPerIntegrationPointValues

static pybind11::object MaterialStateManager_getGradients(
    pybind11::object self) {
  auto& s = self.cast<mgis::behaviour::MaterialStateManager&>();
  return MaterialStateManager_wrapPerIntegrationPointValues(
      s.gradients, s.n, s.gradients_stride, self);
}  // end of MaterialStateManager_getGradients

static pybind11::object MaterialStateManager_getThermodynamicForces(
    pybind11::object self) {
  auto& s = self.cast<mgis::behaviour::MaterialStateManager&>();
  return MaterialStateManager_wrapPerIntegrationPointValues(
      s.thermodynamic_forces, s.n, s.thermodynamic_forces_stride, self);
}  // end of MaterialStateManager_getThermodynamicForces

static pybind11::object MaterialStateManager_getInternalStateVariables(
    pybind11::object self) {
  auto& s = self.cast<mgis::behaviour::MaterialStateManager&>();
  return MaterialStateManager_wrapPerIntegrationPointValues(
      s.internal_state_variables, s.n, s.internal_state_variables_stride, self);
}  // end of MaterialStateManager_getInternalStateVariables

static pybind11::object MaterialStateManager_getStoredEnergies(
    pybind11::object self) {
  auto& s = self.cast<mgis::behaviour::MaterialStateManager&>();
  if (!s.b.computesStoredEnergy) {
    mgis::raise(
        "MaterialStateManager_getStoredEnergies: "
        "the stored energy is not computed by the behaviour");
  }
  return mgis::python::wrapInNumPyArray(s.stored_energies, self);
}  // end of MaterialStateManager_getStoredEnergy

static pybind11::object MaterialStateManager_getDissipatedEnergies(
    pybind11::object self) {
  auto& s = self.cast<mgis::behaviour::MaterialStateManager&>();
  if (!s.b.computesDissipatedEnergy) {
    mgis::raise(
        "MaterialStateManager_getDissipatedEnergies: "
        "the dissipated energy is not computed by the behaviour");
  }
  return mgis::python::wrapInNumPyArray(s.dissipated_energies, self);
}  // end of MaterialStateManager_getDissipatedEnergy

/*!
 * \brief convert values stored externally by a material state manager.
 * The manager keeps the address of the values, so a conversion making a copy
 * is refused. A read-only array is refused unless the values are not updated,
 * since the `update` and `revert` functions write in them.
 * \return the converted array, which must be kept as long as the values are
 * used since it holds the buffer exported by the given object, if any
 * \param[in] o: values
 * \param[in] p: update policy
 */
static pybind11::array_t<double> MaterialStateManager_convertExternalValues(
    const pybind11::object& o,
    const mgis::behaviour::MaterialStateManager::UpdatePolicy p) {
  if (p == mgis::behaviour::MaterialStateManager::NOUPDATE) {
    return mgis::python::mgis_convert_to_array_without_copy(o);
  }
  return mgis::python::mgis_convert_to_mutable_array(o);
}  // end of MaterialStateManager_convertExternalValues

/*!
 * \return the clean-up hook of values stored externally by a material state
 * manager: the hook holds a reference to the given object and releases it.
 * \param[in] o: converted array of the values
 */
static std::function<void()> MaterialStateManager_makeExternalStorageHook(
    pybind11::object o) {
  return [o = std::move(o)]() mutable {
    // the hook may be called by a thread which does not hold the GIL
    const pybind11::gil_scoped_acquire gil;
    o = pybind11::object{};
  };
}  // end of MaterialStateManager_makeExternalStorageHook

static void MaterialStateManager_setMaterialProperty(
    mgis::behaviour::MaterialStateManager& s,
    const std::string& n,
    const mgis::real v) {
  mgis::behaviour::setMaterialProperty(s, n, v);
}  // end of MaterialStateManager_setMaterialProperty

static void MaterialStateManager_setMaterialProperty2(
    mgis::behaviour::MaterialStateManager& s,
    const std::string& n,
    const mgis::real v,
    const mgis::behaviour::MaterialStateManager::UpdatePolicy b) {
  mgis::behaviour::setMaterialProperty(s, n, v, b);
}  // end of MaterialStateManager_setMaterialProperty

static void MaterialStateManager_setMaterialProperty3(
    mgis::behaviour::MaterialStateManager& sm,
    const std::string& n,
    const pybind11::object& o,
    const mgis::behaviour::MaterialStateManager::StorageMode s) {
  if (s == mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE) {
    const auto a = MaterialStateManager_convertExternalValues(
        o, mgis::behaviour::MaterialStateManager::UPDATE);
    setMaterialProperty(sm, n, mgis::python::mgis_convert_to_span(a), s,
                        mgis::behaviour::MaterialStateManager::UPDATE,
                        MaterialStateManager_makeExternalStorageHook(a));
    return;
  }
  setMaterialProperty(sm, n, mgis::python::mgis_convert_to_span(o), s);
}  // end of MaterialStateManager_setMaterialProperty

static void MaterialStateManager_setMaterialProperty4(
    mgis::behaviour::MaterialStateManager& sm,
    const std::string& n,
    const pybind11::object& o,
    const mgis::behaviour::MaterialStateManager::StorageMode s,
    const mgis::behaviour::MaterialStateManager::UpdatePolicy b) {
  if (s == mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE) {
    const auto a = MaterialStateManager_convertExternalValues(o, b);
    setMaterialProperty(sm, n, mgis::python::mgis_convert_to_span(a), s, b,
                        MaterialStateManager_makeExternalStorageHook(a));
    return;
  }
  setMaterialProperty(sm, n, mgis::python::mgis_convert_to_span(o), s, b);
}  // end of MaterialStateManager_setMaterialProperty

static void MaterialStateManager_setMassDensity(
    mgis::behaviour::MaterialStateManager& s, const mgis::real v) {
  mgis::behaviour::setMassDensity(s, v);
}  // end of MaterialStateManager_setMassDensity

static void MaterialStateManager_setMassDensity2(
    mgis::behaviour::MaterialStateManager& s,
    const mgis::real v,
    const mgis::behaviour::MaterialStateManager::UpdatePolicy b) {
  mgis::behaviour::setMassDensity(s, v, b);
}  // end of MaterialStateManager_setMassDensity2

static void MaterialStateManager_setMassDensity3(
    mgis::behaviour::MaterialStateManager& sm,
    const pybind11::object& o,
    const mgis::behaviour::MaterialStateManager::StorageMode s) {
  if (s == mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE) {
    const auto a = MaterialStateManager_convertExternalValues(
        o, mgis::behaviour::MaterialStateManager::UPDATE);
    setMassDensity(sm, mgis::python::mgis_convert_to_span(a), s,
                   mgis::behaviour::MaterialStateManager::UPDATE,
                   MaterialStateManager_makeExternalStorageHook(a));
    return;
  }
  setMassDensity(sm, mgis::python::mgis_convert_to_span(o), s);
}  // end of MaterialStateManager_setMassDensity3

static void MaterialStateManager_setMassDensity4(
    mgis::behaviour::MaterialStateManager& sm,
    const pybind11::object& o,
    const mgis::behaviour::MaterialStateManager::StorageMode s,
    const mgis::behaviour::MaterialStateManager::UpdatePolicy b) {
  if (s == mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE) {
    const auto a = MaterialStateManager_convertExternalValues(o, b);
    setMassDensity(sm, mgis::python::mgis_convert_to_span(a), s, b,
                   MaterialStateManager_makeExternalStorageHook(a));
    return;
  }
  setMassDensity(sm, mgis::python::mgis_convert_to_span(o), s, b);
}  // end of MaterialStateManager_setMassDensity4

static void MaterialStateManager_setExternalStateVariable(
    mgis::behaviour::MaterialStateManager& s,
    const std::string& n,
    const mgis::real v) {
  mgis::behaviour::setExternalStateVariable(s, n, v);
}  // end of MaterialStateManager_setExternalStateVariable

static void MaterialStateManager_setExternalStateVariable2(
    mgis::behaviour::MaterialStateManager& s,
    const std::string& n,
    const mgis::real v,
    const mgis::behaviour::MaterialStateManager::UpdatePolicy b) {
  mgis::behaviour::setExternalStateVariable(s, n, v, b);
}  // end of MaterialStateManager_setExternalStateVariable2

static void MaterialStateManager_setExternalStateVariable3(
    mgis::behaviour::MaterialStateManager& sm,
    const std::string& n,
    const pybind11::object& o,
    const mgis::behaviour::MaterialStateManager::StorageMode s) {
  if (s == mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE) {
    const auto a = MaterialStateManager_convertExternalValues(
        o, mgis::behaviour::MaterialStateManager::UPDATE);
    setExternalStateVariable(sm, n, mgis::python::mgis_convert_to_span(a), s,
                             mgis::behaviour::MaterialStateManager::UPDATE,
                             MaterialStateManager_makeExternalStorageHook(a));
    return;
  }
  setExternalStateVariable(sm, n, mgis::python::mgis_convert_to_span(o), s);
}  // end of MaterialStateManager_setExternalStateVariable3

static void MaterialStateManager_setExternalStateVariable4(
    mgis::behaviour::MaterialStateManager& sm,
    const std::string& n,
    const pybind11::object& o,
    const mgis::behaviour::MaterialStateManager::StorageMode s,
    const mgis::behaviour::MaterialStateManager::UpdatePolicy b) {
  if (s == mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE) {
    const auto a = MaterialStateManager_convertExternalValues(o, b);
    setExternalStateVariable(sm, n, mgis::python::mgis_convert_to_span(a), s, b,
                             MaterialStateManager_makeExternalStorageHook(a));
    return;
  }
  setExternalStateVariable(sm, n, mgis::python::mgis_convert_to_span(o), s, b);
}  // end of MaterialStateManager_setExternalStateVariable4

void declareMaterialStateManager(pybind11::module_&);

void declareMaterialStateManager(pybind11::module_& m) {
  using mgis::size_type;
  using mgis::behaviour::Behaviour;
  using mgis::behaviour::MaterialStateManager;
  using mgis::behaviour::MaterialStateManagerInitializer;
  // documentation of the overloads taking a storage mode
  const auto* const storage_mode_doc =
      "with EXTERNAL_STORAGE, the given array is used without copy: it must "
      "be a contiguous array of float64 values, writable unless the update "
      "policy is NOUPDATE, and the manager keeps a reference to it as long as "
      "it uses it";
  // wrapping the MaterialStateManager::UpdatePolicy enum
  pybind11::enum_<MaterialStateManager::UpdatePolicy>(
      m, "MaterialStateManagerUpdatePolicy")
      .value("UPDATE", MaterialStateManager::UpdatePolicy::UPDATE)
      .value("Update", MaterialStateManager::UpdatePolicy::UPDATE)
      .value("NOUPDATE", MaterialStateManager::UpdatePolicy::NOUPDATE)
      .value("NoUpdate", MaterialStateManager::UpdatePolicy::NOUPDATE);
  // wrapping the MaterialStateManager::StorageMode enum
  pybind11::enum_<MaterialStateManager::StorageMode>(
      m, "MaterialStateManagerStorageMode")
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
  pybind11::class_<MaterialStateManagerInitializer>(
      m, "MaterialStateManagerInitializer")
      .def(pybind11::init<>())
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
  pybind11::class_<MaterialStateManager>(m, "MaterialStateManager")
      .def(pybind11::init<const Behaviour&, const size_type>(),
           pybind11::keep_alive<1, 2>())
      .def(pybind11::init<const Behaviour&, const size_type,
                          const MaterialStateManagerInitializer&>(),
           pybind11::keep_alive<1, 2>(), pybind11::keep_alive<1, 4>())
      .def_readonly("n", &MaterialStateManager::n)
      .def_readonly("number_of_integration_points", &MaterialStateManager::n)
      .def_property_readonly("gradients", &MaterialStateManager_getGradients)
      .def_readonly("gradients_stride", &MaterialStateManager::gradients_stride)
      .def_property_readonly("thermodynamic_forces",
                             &MaterialStateManager_getThermodynamicForces)
      .def_readonly("thermodynamic_forces_stride",
                    &MaterialStateManager::thermodynamic_forces_stride)
      .def_readonly("internal_state_variables_stride",
                    &MaterialStateManager::internal_state_variables_stride)
      .def_property_readonly("stored_energies",
                             &MaterialStateManager_getStoredEnergies)
      .def_property_readonly("dissipated_energies",
                             &MaterialStateManager_getDissipatedEnergies)
      .def_property_readonly("internal_state_variables",
                             &MaterialStateManager_getInternalStateVariables)
      .def("setMaterialProperty", &MaterialStateManager_setMaterialProperty)
      .def("setMaterialProperty", &MaterialStateManager_setMaterialProperty2)
      .def("setMaterialProperty", &MaterialStateManager_setMaterialProperty3,
           storage_mode_doc)
      .def("setMaterialProperty", &MaterialStateManager_setMaterialProperty4,
           storage_mode_doc)
      .def("setMassDensity", &MaterialStateManager_setMassDensity)
      .def("setMassDensity", &MaterialStateManager_setMassDensity2)
      .def("setMassDensity", &MaterialStateManager_setMassDensity3,
           storage_mode_doc)
      .def("setMassDensity", &MaterialStateManager_setMassDensity4,
           storage_mode_doc)
      .def("setExternalStateVariable",
           &MaterialStateManager_setExternalStateVariable)
      .def("setExternalStateVariable",
           &MaterialStateManager_setExternalStateVariable2)
      .def("setExternalStateVariable",
           &MaterialStateManager_setExternalStateVariable3, storage_mode_doc)
      .def("setExternalStateVariable",
           &MaterialStateManager_setExternalStateVariable4, storage_mode_doc);
  // free functions
  m.def("setMaterialProperty", &MaterialStateManager_setMaterialProperty);
  m.def("setMaterialProperty", &MaterialStateManager_setMaterialProperty2);
  m.def("setMaterialProperty", &MaterialStateManager_setMaterialProperty3,
        storage_mode_doc);
  m.def("setMaterialProperty", &MaterialStateManager_setMaterialProperty4,
        storage_mode_doc);
  m.def("setMassDensity", &MaterialStateManager_setMassDensity);
  m.def("setMassDensity", &MaterialStateManager_setMassDensity2);
  m.def("setMassDensity", &MaterialStateManager_setMassDensity3,
        storage_mode_doc);
  m.def("setMassDensity", &MaterialStateManager_setMassDensity4,
        storage_mode_doc);
  m.def("setExternalStateVariable",
        &MaterialStateManager_setExternalStateVariable);
  m.def("setExternalStateVariable",
        &MaterialStateManager_setExternalStateVariable2);
  m.def("setExternalStateVariable",
        &MaterialStateManager_setExternalStateVariable3, storage_mode_doc);
  m.def("setExternalStateVariable",
        &MaterialStateManager_setExternalStateVariable4, storage_mode_doc);

}  // end of declareMaterialStateManager
