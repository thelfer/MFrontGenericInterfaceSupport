/*!
 * \file   bindings/fencis/include/MGIS/FEniCS/Utils.hxx
 * \brief
 * \author Thomas Helfer
 * \date   14/12/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FENICS_UTILS_HXX
#define LIB_MGIS_FENICS_UTILS_HXX

#include "MGIS/FEniCS/Config-FEniCS.hxx"

namespace mgis {

  namespace fenics {

    /*!
     * \return the number of integration points of the finite elements
     * \param[in] e: finite elements
     */
    MGIS_FENICS_EXPORT mgis::size_type getNumberOfIntegrationPoints(
        const dolfin::Mesh&, const dolfin::FiniteElement&);

  }  // end of namespace fenics

}  // end of namespace mgis

#endif /* LIB_MGIS_FENICS_UTILS_HXX */
