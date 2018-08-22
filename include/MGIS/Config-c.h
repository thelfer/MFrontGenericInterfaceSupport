/*!
 * \file   include/MGIS/Config-c.h
 * \brief    
 * \author Thomas Helfer
 * \date   19/06/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_CONFIG_C_H
#define LIB_MGIS_CONFIG_C_H

#ifdef __cplusplus
#include <cstddef>
#else  /*  __cplusplus */
#include <stddef.h>
#endif /*  __cplusplus */

/*!
 * Macro extracted from :
 * "Why is the new C++ visibility support so useful?"
 * from http://gcc.gnu.org/wiki/Visibility
 */
#if defined _WIN32 || defined _WIN64 ||defined __CYGWIN__
#define MGIS_VISIBILITY_IMPORT __declspec(dllimport)
#define MGIS_VISIBILITY_EXPORT __declspec(dllexport)
#define MGIS_VISIBILITY_LOCAL
#else /* defined _WIN32 || defined __CYGWIN__ */
#if (defined __GNUC__) && (! defined __INTEL_COMPILER)
#if __GNUC__ >= 4
#define MGIS_VISIBILITY_IMPORT __attribute__((visibility("default")))
#define MGIS_VISIBILITY_EXPORT __attribute__((visibility("default")))
#define MGIS_VISIBILITY_LOCAL  __attribute__((visibility("hidden")))
#else /* __GNUC__ >= 4 */
#define MGIS_VISIBILITY_IMPORT
#define MGIS_VISIBILITY_EXPORT
#define MGIS_VISIBILITY_LOCAL
#endif /* __GNUC__ >= 4 */
#elif defined __INTEL_COMPILER
#define MGIS_VISIBILITY_IMPORT __attribute__((visibility("default")))
#define MGIS_VISIBILITY_EXPORT __attribute__((visibility("default")))
#define MGIS_VISIBILITY_LOCAL  __attribute__((visibility("hidden")))
#else /* defined __INTEL_COMPILER */
#define MGIS_VISIBILITY_IMPORT
#define MGIS_VISIBILITY_EXPORT
#define MGIS_VISIBILITY_LOCAL
#endif /* defined __INTEL_COMPILER */
#endif /* defined _WIN32 || defined _WIN64 ||defined __CYGWIN__ */

#ifdef MGIS_REAL_TYPE
/*! alias to the numeric type used in the library */
typedef MGIS_REAL_TYPE mgis_real;
#else   /* MGIS_REAL_TYPE */
/*! alias to the numeric type used in the library */
typedef double mgis_real;
#endif /* MGIS_REAL_TYPE */

//! alias to the index type type used in the library
typedef size_t mgis_size_type; 

#define MGIS_INLINE inline
#define MGIS_NORETURN [[noreturn]]

#if defined _WIN32 || defined _WIN64 || defined __CYGWIN__
#if defined MFrontGenericInterface_EXPORTS
#define MGIS_EXPORT MGIS_VISIBILITY_EXPORT
#else
#ifndef MGIS_STATIC_BUILD
#define MGIS_EXPORT MGIS_VISIBILITY_IMPORT
#else
#define MGIS_EXPORT
#endif
#endif
#else
#define MGIS_EXPORT MGIS_VISIBILITY_EXPORT
#endif /* */

#endif /* LIB_MGIS_CONFIG_C_H */
