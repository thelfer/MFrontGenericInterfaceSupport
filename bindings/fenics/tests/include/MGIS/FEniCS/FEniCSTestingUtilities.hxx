/*!
 * \file   bindings/fencis/include/MGIS/FEniCS/FenicsTestingUtilities.hxx
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
 *
 * \note This file contains material strongly inspired by the
 * fenics-solid-materials by Kristian Oelgaard and Garth N. Wells.
 * See <https://bitbucket.org/fenics-apps/fenics-solid-mechanics> for
 * details.
 */

#ifndef LIB_MGIS_FENICS_TESTINGUTILITIES_HXX
#define LIB_MGIS_FENICS_TESTINGUTILITIES_HXX

#include <map>
#include <memory>
#include <dolfin/mesh/SubDomain.h>
#include "MGIS/Config.hxx"

#if defined _WIN32 || defined _WIN64 || defined __CYGWIN__
#if defined MFrontGenericInterfaceFEniCSTestingUtilities_EXPORTS
#define MGIS_FENICS_TESTINGEXPORT MGIS_VISIBILITY_EXPORT
#else
#ifndef MGIS_STATIC_BUILD
#define MGIS_FENICS_TESTINGEXPORT MGIS_VISIBILITY_IMPORT
#else
#define MGIS_FENICS_TESTINGEXPORT
#endif
#endif
#else
#define MGIS_FENICS_TESTINGEXPORT MGIS_VISIBILITY_EXPORT
#endif /* */

namespace mgis {

  namespace fenics {

    /*!
     * \return all the boundaries associated with a unit square
     * The following boundaries are defined:
     * - sx1: defined by the points satisfying x==0
     * - sx2: defined by the points satisfying x==1
     * - sy1: defined by the points satisfying y==0
     * - sy2: defined by the points satisfying y==1
     */
    MGIS_VISIBILITY_EXPORT
    std::map<std::string, std::shared_ptr<dolfin::SubDomain>>
    getUnitSquareBoundaries();
    /*!
     * \return all the boundaries associated with a unit cube
     * The following boundaries are defined:
     * - sx1: defined by the points satisfying x==0
     * - sx2: defined by the points satisfying x==1
     * - sy1: defined by the points satisfying y==0
     * - sy2: defined by the points satisfying y==1
     * - sz1: defined by the points satisfying z==0
     * - sz2: defined by the points satisfying z==1
     */
    MGIS_VISIBILITY_EXPORT
    std::map<std::string, std::shared_ptr<dolfin::SubDomain>>
    getUnitCubeBoundaries();

  }  // end of namespace fenics

}  // end of namespace mgis

#endif /* LIB_MGIS_FENICS_TESTINGUTILITIES_HXX */
