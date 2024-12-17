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
    auto d = MaterialProperty{};
    d.library = l;
    d.material_property = mp;
    d.tfel_version = lm.getTFELVersion(l, mp);
    d.unit_system = lm.getUnitSystem(l, mp);
    d.source = lm.getSource(l, mp);
    d.fct = lm.getMaterialProperty(l, mp);
    d.output = lm.getMaterialPropertyOutputName(l, mp);
    d.inputs = lm.getMaterialPropertyInputsNames(l, mp);
    return d;
  }  // end of load_behaviour

}  // end of namespace mgis::material_property
