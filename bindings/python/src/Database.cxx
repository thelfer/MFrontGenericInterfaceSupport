/*!
 * \file   bindings/python/src/Database.cxx
 * \brief
 * \author Thomas Helfer
 * \date   31/10/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <pybind11/pybind11.h>
#include "MGIS/Database.hxx"

void declareDatabase(pybind11::module_&);

void declareDatabase(pybind11::module_& m) {
  m.def("getDatabase", mgis::getDatabase,
        pybind11::return_value_policy::reference);
}  // end of declareDatabase
