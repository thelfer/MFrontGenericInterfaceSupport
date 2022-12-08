/*!
 * \file   include/MGIS/MaterialProperty/MaterialPropertyFctPtr.hxx
 * \brief
 * \author Thomas Helfer
 * \date   04/10/2022
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_MATERIALPROPERTY_MATERIALPROPERTYFCTPTR_HXX
#define LIB_MGIS_MATERIALPROPERTY_MATERIALPROPERTYFCTPTR_HXX

#ifdef __cplusplus
#include "MGIS/Config.hxx"
#else
#include "MGIS/Config-c.h"
#endif /* __cplusplus */

#include "MGIS/MaterialProperty/OutputStatus.hxx"
#include "MGIS/MaterialProperty/OutOfBoundsPolicy.hxx"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#if defined _WIN32 || defined _WIN64 || defined __CYGWIN__
#define MGIS_MATERIALPROPERTY_ADDCALL_PTR __cdecl*
#else
#define MGIS_MATERIALPROPERTY_ADDCALL_PTR *
#endif

//! \brief a simple alias
typedef mgis_real(
    MGIS_MATERIALPROPERTY_ADDCALL_PTR mgis_mp_MaterialPropertyFctPtr)(
    mgis_mp_OutputStatus* const,       // output status
    const mgis_real* const,            // arguments
    const mgis_size_type,              // number of arguments
    const mgis_mp_OutOfBoundsPolicy);  // out of bounds policy

#ifdef __cplusplus
}  // end of extern "C"
#endif /* __cplusplus */

#ifdef __cplusplus

namespace mgis::material_property {

  //! \brief a simple alias
  using MaterialPropertyFctPtr = ::mgis_mp_MaterialPropertyFctPtr;

}  // end of namespace mgis::material_property

#endif /* __cplusplus */

#endif /* LIB_MGIS_MATERIALPROPERTY_MATERIALPROPERTYFCTPTR_HXX */
