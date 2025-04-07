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

#include <pybind11/pybind11.h>

// forward declaration
void declareModel(pybind11::module_&);

PYBIND11_MODULE(model, m) { declareModel(m); }  // end of module model
