/*!
 * \file   VectorConverter.hxx
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

#ifndef LIB_MGIS_PYTHON_VECTORCONVERTER_HXX
#define LIB_MGIS_PYTHON_VECTORCONVERTER_HXX

#include <iostream>
#include <vector>

#ifdef _MSC_VER
#ifdef and
#undef and
#endif
#ifdef xor
#undef xor
#endif
#ifdef or
#undef or
#endif
#endif /* _MSC_VER */

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

namespace mgis {

  namespace python {

    template <typename T>
    boost::python::list convert_vector_to_list(const std::vector<T>& v) {
      boost::python::list l;
      for (const auto& e : v) {
        l.append(e);
      }
      return l;
    }  // end of convert

    /*!
     * convert vector to python list
     */
    template <typename T>
    struct vector_to_python_list {
      static PyObject* convert(const T& o) {
        boost::python::list l;
        for (const auto& v : o) {
          l.append(v);
        }
        return boost::python::incref(l.ptr());
      }
    };

    template <typename T>
    struct vector_from_python_list {
      vector_from_python_list() {
        using boost::python::type_id;
        using namespace boost::python::converter;
        registry::push_back(&convertible, &construct, type_id<T>());
      }

      static void* convertible(PyObject* ptr) {
        using namespace boost::python;
        if (!PyList_Check(ptr)) {
          return nullptr;
        }
        handle<> h(borrowed(ptr));
        list l(h);
        stl_input_iterator<object> p(l);
        stl_input_iterator<object> pe;
        while (p != pe) {
          extract<typename T::value_type> e(*p);
          if (!e.check()) {
            return nullptr;
          }
          ++p;
        }
        return ptr;
      }

      static void construct(
          PyObject* ptr,
          boost::python::converter::rvalue_from_python_stage1_data* data) {
        using namespace boost::python;
        using namespace boost::python::converter;
        typedef rvalue_from_python_storage<T> py_storage;
        assert(PyList_Check(ptr));
        handle<> h(borrowed(ptr));
        list l(h);
        stl_input_iterator<object> p(l);
        stl_input_iterator<object> pe;
        // Grab pointer to memory into which to construct the new T
        void* storage = reinterpret_cast<py_storage*>(data)->storage.bytes;
        // in-place construct the new T using the character data
        // extraced from the python object
        new (storage) T();
        T& v = *(static_cast<T*>(storage));
        while (p != pe) {
          extract<typename T::value_type> e(*p);
          assert(e.check());
          v.push_back(e());
          ++p;
        }
        // Stash the memory chunk pointer for later use by boost.python
        data->convertible = storage;
      }
    };

    template <typename T>
    static void initializeVectorConverter() {
      using namespace boost::python;
      // register the to-python converter
      to_python_converter<T, vector_to_python_list<T>>();
      // register the from-python converter
      vector_from_python_list<T>();
    }

  }  // end of namespace python

}  // end of namespace mgis

#endif /* LIB_MGIS_PYTHON_VECTORCONVERTER_HXX */
