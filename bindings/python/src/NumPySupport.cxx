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
#include "MGIS/Raise.hxx"
#include "MGIS/Python/NumPySupport.hxx"

namespace mgis {

  namespace python {

#if PY_MAJOR_VERSION == 2
    static void wrapInitializeNumPy() { import_array(); }
#else
    static void* wrapInitializeNumPy() { import_array(); return nullptr;}
#endif

    void initializeNumPy() { wrapInitializeNumPy(); }  // end of initializeNumPy

    boost::python::object wrapInNumPyArray(mgis::span<double>& v) {
      npy_intp dims[1] = {v.size()};
      auto* const arr =
          PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, v.data());
      boost::python::handle<> handle(arr);
      return boost::python::object(handle);
    }  // end of wrapInNumPyArray

    boost::python::object wrapInNumPyArray(std::vector<double>& v) {
      auto view = mgis::span<double>(v);
      return wrapInNumPyArray(view);
    }  // end of wrapInNumPyArray

    boost::python::object wrapInNumPyArray(
        mgis::variant<mgis::span<double>, std::vector<double>>& v) {
      if (mgis::holds_alternative<std::vector<double>>(v)) {
        return wrapInNumPyArray(mgis::get<std::vector<double>>(v));
      }
      return wrapInNumPyArray(mgis::get<mgis::span<double>>(v));
    }  // end of wrapInNumPyArray

    boost::python::object wrapInNumPyArray(mgis::span<double>& v,
                                           const size_type nc) {
      npy_intp dims[2] = {static_cast<npy_intp>(v.size() / nc),
			  static_cast<npy_intp>(nc)};
      auto* const arr =
          PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, v.data());
      boost::python::handle<> handle(arr);
      return boost::python::object(handle);
    }  // end of wrapInNumPyArray

    boost::python::object wrapInNumPyArray(std::vector<double>& v,
                                           const size_type nc) {
      auto view = mgis::span<double>(v);
      return wrapInNumPyArray(view, nc);
    }  // end of wrapInNumPyArray

    boost::python::object wrapInNumPyArray(
        mgis::variant<mgis::span<double>, std::vector<double>>& v,
        const size_type nc) {
      if (mgis::holds_alternative<std::vector<double>>(v)) {
        return wrapInNumPyArray(mgis::get<std::vector<double>>(v), nc);
      }
      return wrapInNumPyArray(mgis::get<mgis::span<double>>(v), nc);
    }  // end of wrapInNumPyArray

    boost::python::object wrapInNumPyArray(mgis::span<double>& v,
                                           const size_type nl,
                                           const size_type nc) {
      npy_intp dims[3] = {static_cast<npy_intp>(v.size() / (nc * nl)),
			  static_cast<npy_intp>(nl),
			  static_cast<npy_intp>(nc)};
      auto* const arr =
          PyArray_SimpleNewFromData(3, dims, NPY_DOUBLE, v.data());
      boost::python::handle<> handle(arr);
      return boost::python::object(handle);
    }  // end of wrapInNumPyArray

    boost::python::object wrapInNumPyArray(std::vector<double>& v,
                                           const size_type nl,
                                           const size_type nc) {
      auto view = mgis::span<double>(v);
      return wrapInNumPyArray(view, nl, nc);
    }  // end of wrapInNumPyArray

    boost::python::object wrapInNumPyArray(
        mgis::variant<mgis::span<double>, std::vector<double>>& v,
        const size_type nl,
        const size_type nc) {
      if (mgis::holds_alternative<std::vector<double>>(v)) {
        return wrapInNumPyArray(mgis::get<std::vector<double>>(v), nl, nc);
      }
      return wrapInNumPyArray(mgis::get<mgis::span<double>>(v), nl, nc);
    }  // end of wrapInNumPyArray

    mgis::span<mgis::real> mgis_convert_to_span(
        const boost::python::object& o) {
      PyArrayObject* a = nullptr;
      if (!PyArg_ParseTuple(o.ptr(), "O!", &PyArray_Type, &a)) {
        mgis::raise(
            "convert_to_span: argument not convertible to PyArrayObject");
      }
      if (a->descr->type_num != NPY_DOUBLE) {
        mgis::raise("convert_to_span: invalid numpy object");
      }
      if (a == nullptr) {
        mgis::raise(
            "convert_to_span: argument not convertible to PyArrayObject");
      }
      // number of dimensions
      const auto nd = PyArray_NDIM(a);
      if (nd != 1) {
        mgis::raise("convert_to_span: expected one dimensional array");
      }
      auto* const shape = PyArray_DIMS(a);
      const auto n = shape[0];
      auto* const values = static_cast<double*>(PyArray_GETPTR1(a, 0));
      return {values, static_cast<mgis::span<mgis::real>::index_type>(n)};
    }  // end of mgis_convert_to_span

  }  // end of namespace python

}  // end of namespace mgis
