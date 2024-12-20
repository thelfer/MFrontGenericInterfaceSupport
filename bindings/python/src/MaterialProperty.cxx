/*!
 * \file   bindings/python/src/MaterialProperty.cxx
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
#include "MGIS/MaterialProperty/MaterialProperty.hxx"

void declareMaterialProperty(pybind11::module_&);

void declareMaterialProperty(pybind11::module_& m) {
  using mgis::material_property::MaterialProperty;
  // wrapping the MaterialProperty class
  pybind11::class_<MaterialProperty>(m, "MaterialProperty")
      .def_readonly(
          "library", &MaterialProperty::library,
          "name of the library in which the material property is implemented")
      .def_readonly("material_property", &MaterialProperty::material_property,
                    "name of the material property")
      .def_readonly("source", &MaterialProperty::source,
                    "name of the `MFront` source file")
      .def_readonly("tfel_version", &MaterialProperty::tfel_version,
                    "version of TFEL used to generate the material property")
      .def_readonly("output", &MaterialProperty::output,
                    "output of the material property")
      .def_readonly("inputs", &MaterialProperty::inputs,
                    "inputs of the material property");
  // wrapping free functions
  m.def("load", mgis::material_property::load);

}  // end of declareMaterialProperty
