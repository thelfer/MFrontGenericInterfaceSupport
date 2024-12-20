/*!
 * \file   bindings/python/src/Model.cxx
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
#include "MGIS/Model/Model.hxx"

// forward declaration
void declareModel(pybind11::module_& m);

void declareModel(pybind11::module_& m) {
  // wrapping free functions
  m.def("load", mgis::model::load);
}  // end of declareModel
