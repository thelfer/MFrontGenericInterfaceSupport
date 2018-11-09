/*!
 * \file   NumPySupport.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   08/11/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <numpy/ndarrayobject.h>
#include "MGIS/Python/NumPySupport.hxx"

namespace mgis {

  namespace python {

#if PY_MAJOR_VERSION == 2
    static void wrapInitializeNumPy() { import_array(); }
#else
    static void* wrapInitializeNumPy() { import_array(); }
#endif

    void initializeNumPy() { wrapInitializeNumPy(); }  // end of initializeNumPy

    boost::python::object wrapInNumPyArray(std::vector<double>& v) {
      npy_intp dims[1] = {v.size()};
      auto* const arr =
          PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, v.data());
      boost::python::handle<> handle(arr);
      return boost::python::object(handle);
    }  // end of wrapInNumPyArray

  }  // end of namespace python

}  // end of namespace mgis
