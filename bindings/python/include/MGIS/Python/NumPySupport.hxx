/*!
 * \file   bindings/python/include/MGIS/Python/NumPySupport.hxx
 * \brief
 * \author Thomas Helfer
 * \date   07/11/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_PYTHON_NUMPYSUPPORT_HXX
#define LIB_MGIS_PYTHON_NUMPYSUPPORT_HXX

#include <numpy/ndarrayobject.h>
#include <vector>
#include <boost/python/object.hpp>

namespace mgis {

  namespace python {

#if PY_MAJOR_VERSION == 2
    /*!
     * \brief initialize NumPy
     */
    inline void initializeNumPy() { import_array(); }
#else
    inline void* initializeNumPy() { import_array(); }
#endif

    /*!
     * \brief create a ndarray object from a vector.
     * The ndarray does not own the data, the lifetime of which is handled by
     * the vector.
     * \param[in] v: vector holding the values
     */
    inline boost::python::object wrapInNumPyArray(std::vector<double>& v) {
      npy_intp dims[1] = {v.size()};
      auto* const arr = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, &v[0]);
      boost::python::handle<> handle(arr);
      return boost::python::object(handle);
    }  // end of wrapInNumPyArray

  }  // end of namespace python

}  // end of namespace mgis

#endif /* LIB_MGIS_PYTHON_NUMPYSUPPORT_HXX */
