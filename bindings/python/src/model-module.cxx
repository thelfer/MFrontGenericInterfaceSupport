/*!
 * \file   bindings/python/src/model-module.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/11/2021
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <string>
#include <pybind11/pybind11.h>

// forward declaration
void declareModel(pybind11::module_&);

PYBIND11_MODULE(model, m) {
  // the types of the signatures of this module are declared by the behaviour
  // module, which shall thus be loaded first, so that they appear by their
  // python names. The name of the package depends on the suffix of MGIS.
  const auto name = m.attr("__name__").cast<std::string>();
  pybind11::module_::import(
      (name.substr(0, name.rfind('.')) + ".behaviour").c_str());
  declareModel(m);
}  // end of module model
