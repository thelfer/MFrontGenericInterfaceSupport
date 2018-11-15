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

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include "MGIS/Python/NumPySupport.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"

void declareMaterialDataManager();

static boost::python::object MaterialDataManager_getK(
    mgis::behaviour::MaterialDataManager& d) {
  const auto nl = d.s0.gradients_stride;
  const auto nc = d.s1.thermodynamic_forces_stride;
  return mgis::python::wrapInNumPyArray(d.K, nl, nc);
}  // end of MaterialDataManager_getK

void declareMaterialDataManager() {
  using mgis::size_type;
  using mgis::behaviour::Behaviour;
  using mgis::behaviour::MaterialDataManager;
  // pointers to free functions to disambiguate the function resolution
  void (*ptr_update)(MaterialDataManager&) = &mgis::behaviour::update;
  void (*ptr_revert)(MaterialDataManager&) = &mgis::behaviour::revert;
  // exporting the BehaviourData class
  boost::python::class_<MaterialDataManager>(
      "MaterialDataManager",
      boost::python::init<const Behaviour&, const size_type>())
      .def_readonly("n", &MaterialDataManager::n)
      .def_readonly("number_of_integration_points", &MaterialDataManager::n)
      .add_property("s0", &MaterialDataManager::s0)
      .add_property("s1", &MaterialDataManager::s1)
      .add_property("K", &MaterialDataManager_getK)
      .def("update", ptr_update)
      .def("revert", ptr_revert);
  // free functions
  boost::python::def("update", ptr_update);
  boost::python::def("revert", ptr_revert);

}  // end of declareMaterialDataManager