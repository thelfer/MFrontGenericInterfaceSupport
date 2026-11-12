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

#include <span>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "MGIS/Config.hxx"

namespace mgis::python {

  /*!
   * \brief create a 1D-ndarray object from a vector.
   * The ndarray does not own the data, the lifetime of which is handled by
   * the vector.
   * \param[in] v: values
   * \param[in] base: object owning the values, kept alive by the ndarray
   */
  pybind11::array_t<double> wrapInNumPyArray(std::span<double>&,
                                             pybind11::handle);
  /*!
   * \brief create a 1D-ndarray object from a vector.
   * The ndarray does not own the data, the lifetime of which is handled by
   * the vector.
   * \param[in] v: vector holding the values
   * \param[in] base: object owning the values, kept alive by the ndarray
   */
  pybind11::array_t<double> wrapInNumPyArray(std::vector<double>&,
                                             pybind11::handle);
  /*!
   * \brief create a 2D-ndarray object from a vector.
   * The ndarray does not own the data, the lifetime of which is handled by
   * the vector.
   * \param[in] v: values
   * \param[in] nc: number of columns
   * \param[in] base: object owning the values, kept alive by the ndarray
   */
  pybind11::array_t<double> wrapInNumPyArray(std::span<double>&,
                                             const mgis::size_type,
                                             pybind11::handle);
  /*!
   * \brief create a 2D-ndarray object from a vector.
   * The ndarray does not own the data, the lifetime of which is handled by
   * the vector.
   * \param[in] v: vector holding the values
   * \param[in] nc: number of columns
   * \param[in] base: object owning the values, kept alive by the ndarray
   */
  pybind11::array_t<double> wrapInNumPyArray(std::vector<double>&,
                                             const mgis::size_type,
                                             pybind11::handle);
  /*!
   * \brief create a 3D-ndarray object from a vector.
   * The ndarray does not own the data, the lifetime of which is handled by
   * the vector.
   * \param[in] v: values
   * \param[in] nl: number of line
   * \param[in] nc: number of columns
   * \param[in] base: object owning the values, kept alive by the ndarray
   */
  pybind11::array_t<double> wrapInNumPyArray(std::span<double>&,
                                             const mgis::size_type,
                                             const mgis::size_type,
                                             pybind11::handle);
  /*!
   * \brief create a 3D-ndarray object from a vector.
   * The ndarray does not own the data, the lifetime of which is handled by
   * the vector.
   * \param[in] v: vector holding the values
   * \param[in] nl: number of line
   * \param[in] nc: number of columns
   * \param[in] base: object owning the values, kept alive by the ndarray
   */
  pybind11::array_t<double> wrapInNumPyArray(std::vector<double>&,
                                             const mgis::size_type,
                                             const mgis::size_type,
                                             pybind11::handle);

  std::span<mgis::real> mgis_convert_to_span(const pybind11::array_t<double>&);
  /*!
   * \brief convert an object into an array holding the buffer exported by this
   * object, if any. An exception is thrown if the object is converted into a
   * copy.
   * \param[in] o: object
   */
  pybind11::array_t<double> mgis_convert_to_array_without_copy(
      const pybind11::object&);
  /*!
   * \brief convert an object in which values are written into a span.
   * Contrary to `mgis_convert_to_span`, an exception is thrown if the object
   * is converted into a copy, in which the written values would be lost, or
   * if the object is read-only.
   * \param[in] o: object
   */
  std::span<mgis::real> mgis_convert_to_mutable_span(const pybind11::object&);
  /*!
   * \brief convert an object in which values are written into an array holding
   * the buffer exported by this object, if any. Contrary to the object, this
   * array protects the values from a reallocation.
   * An exception is thrown if the object is converted into a copy, or if the
   * object is read-only.
   * \param[in] o: object
   */
  pybind11::array_t<double> mgis_convert_to_mutable_array(
      const pybind11::object&);
  /*!
   * \brief convert an object in which values are written into a span used
   * after the call: the converted array is kept alive by the given object.
   * \param[in] o: object
   * \param[in] u: object using the span
   */
  std::span<mgis::real> mgis_convert_to_mutable_span(const pybind11::object&,
                                                     pybind11::handle);
  /*!
   * \brief convert an object given to an initializer: same as
   * `mgis_convert_to_mutable_span`, but an exception is thrown if the object is
   * empty, as an empty array means that no array is given to the initializer,
   * so the object would be silently ignored.
   * \param[in] o: object
   * \param[in] i: initializer
   */
  std::span<mgis::real> mgis_convert_to_bound_span(const pybind11::object&,
                                                   pybind11::handle);

}  // end of namespace mgis::python

#endif /* LIB_MGIS_PYTHON_NUMPYSUPPORT_HXX */
