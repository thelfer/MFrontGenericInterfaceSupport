/*!
 * \file   ThreadPool.cxx
 * \brief
 * \author Thomas Helfer
 * \date   12/11/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <pybind11/pybind11.h>
#include "MGIS/ThreadPool.hxx"

void declareThreadPool(pybind11::module_& m) {
  pybind11::class_<mgis::ThreadPool>(m, "ThreadPool")
      .def(pybind11::init<mgis::size_type>());
}  // end of declareThreadPool
