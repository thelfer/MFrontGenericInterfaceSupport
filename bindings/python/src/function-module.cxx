/*!
 * \file   bindings/python/src/function-module.cxx
 * \brief
 * \author Thomas Helfer
 * \date   18/11/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <pybind11/pybind11.h>

// forward declarations
void declareFunction(pybind11::module_&);

PYBIND11_MODULE(function, m) { declareFunction(m); }  // end of module function
