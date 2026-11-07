/*!
 * \file   tests/GenericParallelInvalidMaterialProperties.cxx
 * \brief  This file exports the symbols of material properties of the
 * `generic-parallel` interface which `MGIS` must refuse: functions of a
 * version which is not supported, invalid memory space, library generated
 * by a version of the interface which does not export the version of its
 * functions. No function is defined: those material properties are refused
 * before any function is loaded.
 * \date   24/10/2026
 * \copyright Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include "MGIS/Config-c.h"
#include "MGIS/MaterialProperty/BatchMaterialPropertyFctPtr.hxx"

extern "C" {

//! \brief version of the functions which is not supported
MGIS_VISIBILITY_EXPORT const char* UnsupportedVersion_mfront_interface =
    "GenericParallel";
MGIS_VISIBILITY_EXPORT unsigned short
    UnsupportedVersion_generic_parallel_api_version =
        MGIS_MP_GENERIC_PARALLEL_API_VERSION + 1;
MGIS_VISIBILITY_EXPORT const char* UnsupportedVersion_generic_parallel_backend =
    "stlpar";
MGIS_VISIBILITY_EXPORT const char*
    UnsupportedVersion_generic_parallel_memory_space = "host";

//! \brief invalid memory space
MGIS_VISIBILITY_EXPORT const char* InvalidMemorySpace_mfront_interface =
    "GenericParallel";
MGIS_VISIBILITY_EXPORT unsigned short
    InvalidMemorySpace_generic_parallel_api_version =
        MGIS_MP_GENERIC_PARALLEL_API_VERSION;
MGIS_VISIBILITY_EXPORT const char* InvalidMemorySpace_generic_parallel_backend =
    "stlpar";
MGIS_VISIBILITY_EXPORT const char*
    InvalidMemorySpace_generic_parallel_memory_space = "nowhere";

//! \brief no symbol giving the version of the functions
MGIS_VISIBILITY_EXPORT const char* NoVersion_mfront_interface =
    "GenericParallel";

}  // end of extern "C"
