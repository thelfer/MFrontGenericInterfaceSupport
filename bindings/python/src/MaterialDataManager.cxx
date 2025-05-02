/*!
 * \file   MaterialDataManager.cxx
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
#include "MGIS/Python/NumPySupport.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"

void declareMaterialDataManager(pybind11::module_&);

static void MaterialDataManagerInitializer_bindTangentOperator(
    mgis::behaviour::MaterialDataManagerInitializer& i, pybind11::object K) {
  i.K = mgis::python::mgis_convert_to_span(K);
}  // end of MaterialDataManagerInitializer_bindTangentOperator

static void MaterialDataManagerInitializer_bindSpeedOfSound(
    mgis::behaviour::MaterialDataManagerInitializer& i, pybind11::object vs) {
  i.speed_of_sound = mgis::python::mgis_convert_to_span(vs);
}  // end of MaterialDataManagerInitializer_bindSpeedOfSound

static void MaterialDataManager_useExternalArrayOfTangentOperatorBlocks(
    mgis::behaviour::MaterialDataManager& m, pybind11::object K) {
  m.useExternalArrayOfTangentOperatorBlocks(
      mgis::python::mgis_convert_to_span(K));
}  // end of MaterialDataManager_useExternalArrayOfTangentOperatorBlocks

static void MaterialDataManager_useExternalArrayOfSpeedOfSounds(
    mgis::behaviour::MaterialDataManager& m, pybind11::object vs) {
  m.useExternalArrayOfSpeedOfSounds(mgis::python::mgis_convert_to_span(vs));
}  // end of MaterialDataManager_useExternalArrayOfSpeedOfSounds

static pybind11::object MaterialDataManager_getK(
    mgis::behaviour::MaterialDataManager& d) {
  if (d.b.to_blocks.size() == 1u) {
    const auto nl =
        getVariableSize(d.b.to_blocks.front().second, d.b.hypothesis);
    const auto nc =
        getVariableSize(d.b.to_blocks.front().first, d.b.hypothesis);
    return mgis::python::wrapInNumPyArray(d.K, nc, nl);
  }
  const auto s = getTangentOperatorArraySize(d.b);
  return mgis::python::wrapInNumPyArray(d.K, s);
}  // end of MaterialDataManager_getK

void declareMaterialDataManager(pybind11::module_& m) {
  using mgis::size_type;
  using mgis::behaviour::Behaviour;
  using mgis::behaviour::MaterialDataManager;
  using mgis::behaviour::MaterialDataManagerInitializer;
  // pointers to free functions to disambiguate the function resolution
  void (*ptr_update)(MaterialDataManager&) = &mgis::behaviour::update;
  void (*ptr_revert)(MaterialDataManager&) = &mgis::behaviour::revert;
  // exporting the MaterialDataManager class
  pybind11::class_<MaterialDataManagerInitializer>(
      m, "MaterialDataManagerInitializer")
      .def_readonly("s0", &MaterialDataManagerInitializer::s0)
      .def_readonly("s1", &MaterialDataManagerInitializer::s1)
      .def("bindTangentOperator",
           &MaterialDataManagerInitializer_bindTangentOperator,
           "use the given array to store the tangent operator blocks")
      .def("bindSpeedOfSound", &MaterialDataManagerInitializer_bindSpeedOfSound,
           "use the given array to store the speed of sounds");
  // exporting the MaterialDataManager class
  pybind11::class_<MaterialDataManager>(m, "MaterialDataManager")
      .def(pybind11::init<const Behaviour&, const size_type>())
      .def(pybind11::init<const Behaviour&, const size_type,
                          const MaterialDataManagerInitializer&>())
      .def("setThreadSafe", &MaterialDataManager::setThreadSafe,
           "specify if various operations performed by the MaterialDataManager "
           "(memory allocations) shall be protected by a mutex")
      .def("allocateArrayOfTangentOperatorBlocks",
           &MaterialDataManager::allocateArrayOfTangentOperatorBlocks,
           "allocate the tangent operator blocks")
      .def("useExternalArrayOfTangentOperatorBlocks",
           &MaterialDataManager_useExternalArrayOfTangentOperatorBlocks,
           "use and externally defined array to store the tangent operator "
           "blocks")
      .def("releaseArrayOfTangentOperatorBlocks",
           &MaterialDataManager::releaseArrayOfTangentOperatorBlocks,
           "release the arrays of tangent operator blocks")
      .def("allocateArrayOfSpeedOfSounds",
           &MaterialDataManager::allocateArrayOfSpeedOfSounds,
           "allocate the array of speed of sounds")
      .def("useExternalArrayOfSpeedOfSounds",
           &MaterialDataManager_useExternalArrayOfSpeedOfSounds,
           "use and externally defined array to store the spped of sounds")
      .def("releaseArrayOfSpeedOfSounds",
           &MaterialDataManager::releaseArrayOfSpeedOfSounds,
           "release the array of speed of sounds")
      .def_readonly("n", &MaterialDataManager::n)
      .def_readonly("number_of_integration_points", &MaterialDataManager::n)
      .def_readonly("s0", &MaterialDataManager::s0)
      .def_readonly("s1", &MaterialDataManager::s1)
      .def_property_readonly("K", &MaterialDataManager_getK)
      .def("update", ptr_update)
      .def("revert", ptr_revert);
  // free functions
  m.def("update", ptr_update);
  m.def("revert", ptr_revert);

}  // end of declareMaterialDataManager
