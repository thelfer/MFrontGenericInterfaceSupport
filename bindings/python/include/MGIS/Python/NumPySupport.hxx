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

namespace mgis {

  namespace python {

    //! \brief initialize NumPy
    void initializeNumPy();

    /*!
     * \brief create a ndarray object from a vector.
     * The ndarray does not own the data, the lifetime of which is handled by
     * the vector.
     * \param[in] v: vector holding the values
     */
    boost::python::object wrapInNumPyArray(std::vector<double>&);

  }  // end of namespace python

}  // end of namespace mgis

#endif /* LIB_MGIS_PYTHON_NUMPYSUPPORT_HXX */
