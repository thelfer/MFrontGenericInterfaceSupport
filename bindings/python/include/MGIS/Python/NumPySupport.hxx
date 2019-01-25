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

#include <vector>
#include <boost/python/object.hpp>
#include "MGIS/Config.hxx"
#include "MGIS/Variant.hxx"
#include "MGIS/Span.hxx"

namespace mgis {

  namespace python {

    //! \brief initialize NumPy
    void initializeNumPy();

    /*!
     * \brief create a 1D-ndarray object from a vector.
     * The ndarray does not own the data, the lifetime of which is handled by
     * the vector.
     * \param[in] v: values
     */
    boost::python::object wrapInNumPyArray(mgis::span<double>&);
    /*!
     * \brief create a 1D-ndarray object from a vector.
     * The ndarray does not own the data, the lifetime of which is handled by
     * the vector.
     * \param[in] v: values
     */
    boost::python::object wrapInNumPyArray(
        mgis::variant<mgis::span<double>, std::vector<double>>&);
    /*!
     * \brief create a 1D-ndarray object from a vector.
     * The ndarray does not own the data, the lifetime of which is handled by
     * the vector.
     * \param[in] v: vector holding the values
     */
    boost::python::object wrapInNumPyArray(std::vector<double>&);
    /*!
     * \brief create a 2D-ndarray object from a vector.
     * The ndarray does not own the data, the lifetime of which is handled by
     * the vector.
     * \param[in] v: values
     * \param[in] nc: number of columns
     */
    boost::python::object wrapInNumPyArray(mgis::span<double>&,
                                           const mgis::size_type);
    /*!
     * \brief create a 2D-ndarray object from a vector.
     * The ndarray does not own the data, the lifetime of which is handled by
     * the vector.
     * \param[in] v: values
     * \param[in] nc: number of columns
     */
    boost::python::object wrapInNumPyArray(
        mgis::variant<mgis::span<double>, std::vector<double>>&,
        const mgis::size_type);
    /*!
     * \brief create a 2D-ndarray object from a vector.
     * The ndarray does not own the data, the lifetime of which is handled by
     * the vector.
     * \param[in] v: vector holding the values
     * \param[in] nc: number of columns
     */
    boost::python::object wrapInNumPyArray(std::vector<double>&,
                                           const mgis::size_type);
    /*!
     * \brief create a 3D-ndarray object from a vector.
     * The ndarray does not own the data, the lifetime of which is handled by
     * the vector.
     * \param[in] v: values
     * \param[in] nl: number of line
     * \param[in] nc: number of columns
     */
    boost::python::object wrapInNumPyArray(mgis::span<double>&,
                                           const mgis::size_type,
                                           const mgis::size_type);
    /*!
     * \brief create a 3D-ndarray object from a vector.
     * The ndarray does not own the data, the lifetime of which is handled by
     * the vector.
     * \param[in] v: values
     * \param[in] nl: number of line
     * \param[in] nc: number of columns
     */
    boost::python::object wrapInNumPyArray(
        mgis::variant<mgis::span<double>, std::vector<double>>&,
        const mgis::size_type,
        const mgis::size_type);
    /*!
     * \brief create a 3D-ndarray object from a vector.
     * The ndarray does not own the data, the lifetime of which is handled by
     * the vector.
     * \param[in] v: vector holding the values
     * \param[in] nl: number of line
     * \param[in] nc: number of columns
     */
    boost::python::object wrapInNumPyArray(std::vector<double>&,
                                           const mgis::size_type,
                                           const mgis::size_type);

    mgis::span<mgis::real> mgis_convert_to_span(const boost::python::object&);

  }  // end of namespace python

}  // end of namespace mgis

#endif /* LIB_MGIS_PYTHON_NUMPYSUPPORT_HXX */
