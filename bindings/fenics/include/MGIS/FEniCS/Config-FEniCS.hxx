/*!
 * \file   bindings/fenics/include/MGIS/FEniCS/Config-FEniCS.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   16/12/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FENICS_CONFIG_HXX
#define LIB_MGIS_FENICS_CONFIG_HXX

#include "MGIS/Config.hxx"

#if defined _WIN32 || defined _WIN64 || defined __CYGWIN__
#if defined MFrontGenericInterface_FEniCS_EXPORTS
#define MGIS_FENICS_EXPORT MGIS_VISIBILITY_EXPORT
#else
#ifndef MGIS_STATIC_BUILD
#define MGIS_FENICS_EXPORT MGIS_VISIBILITY_IMPORT
#else
#define MGIS_FENICS_EXPORT
#endif
#endif
#else
#define MGIS_FENICS_EXPORT MGIS_VISIBILITY_EXPORT
#endif /* */

// forward declarations
namespace dolfin {
  class Cell;
  class FiniteElement;
  class Function;
  class Mesh;
} // end of namespace dolphin

#endif /* LIB_MGIS_FENICS_CONFIG_HXX */
