/*!
 * \file   src/MaterialProperty.cxx
 * \brief
 * \author Thomas Helfer
 * \date   04/10/2022
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <string>
#include "MGIS/Raise.hxx"
#include "MGIS/LibrariesManager.hxx"
#include "MGIS/MaterialProperty/MaterialProperty.hxx"

namespace mgis::material_property {

  MaterialProperty::MaterialProperty() = default;
  MaterialProperty::MaterialProperty(MaterialProperty &&) = default;
  MaterialProperty::MaterialProperty(const MaterialProperty &) = default;
  MaterialProperty &MaterialProperty::operator=(MaterialProperty &&) = default;
  MaterialProperty &MaterialProperty::operator=(const MaterialProperty &) =
      default;
  MaterialProperty::~MaterialProperty() = default;

  MaterialProperty load(const std::string &l, const std::string &mp) {
    auto &lm = mgis::LibrariesManager::get();
    auto raise_if = [&l, &mp](const bool c, const std::string &msg) {
      if (c) {
        mgis::raise("load: " + msg +
                    ".\nError while trying to load material property '" + mp +
                    "' in library '" + l + "'\n");
      }
    };
    auto d = MaterialProperty{};
    d.library = l;
    d.material_property = mp;
    d.tfel_version = lm.getTFELVersion(l, mp);
    d.unit_system = lm.getUnitSystem(l, mp);
    d.source = lm.getSource(l, mp);
    d.mfront_interface = lm.getInterface(l, mp);
    if (lm.hasGenericParallelFunctions(l, mp)) {
      // functions of batches of points (`generic-parallel` interface)
      d.batch_api_version = lm.getGenericParallelAPIVersion(l, mp);
      raise_if(d.batch_api_version != MGIS_MP_GENERIC_PARALLEL_API_VERSION,
               "unmatched version of the functions of the generic-parallel "
               "interface (the material property uses version " +
                   std::to_string(d.batch_api_version) +
                   ", mgis uses version " +
                   std::to_string(MGIS_MP_GENERIC_PARALLEL_API_VERSION) + ")");
      d.backend = lm.getGenericParallelBackend(l, mp);
      d.memory_space = lm.getGenericParallelMemorySpace(l, mp);
      raise_if((d.memory_space != "host") && (d.memory_space != "device") &&
                   (d.memory_space != "any"),
               "invalid memory space '" + d.memory_space + "'");
      d.execution_policy = lm.getGenericParallelExecutionPolicy(l, mp);
      d.batch_fct = lm.getBatchMaterialProperty(l, mp);
      d.batch_report_fct = lm.getBatchMaterialPropertyReportFunction(l, mp);
    } else {
      d.fct = lm.getMaterialProperty(l, mp);
    }
    d.output = lm.getMaterialPropertyOutputName(l, mp);
    d.inputs = lm.getMaterialPropertyInputsNames(l, mp);
    return d;
  }  // end of load

  bool isGenericParallelMaterialProperty(const MaterialProperty &mp) noexcept {
    return mp.batch_fct != nullptr;
  }  // end of isGenericParallelMaterialProperty

}  // end of namespace mgis::material_property
