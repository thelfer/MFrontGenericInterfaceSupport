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
#else /*  __cplusplus */
#include <stddef.h>
#endif /*  __cplusplus */

/*!
 * \brief define the API used by the behaviour support.
 *
 * The API mostly refers to the various functions provided by `MGIS` but also
 * the data strutures used to manage behaviours (See `BehaviourData`,
 * `MaterialDataManger`, etc...).
 *
 * - 0 corresponds to the initial api provided by the TFEL project
 *   (versions 3.2.x, versions 3.3.x and 3.4.0)`.
 * - 1 corresponds the initial api provided by the TFEL project (versions 3.4.1
 *   and after that and 4.x).
 *
 */
#define MGIS_BEHAVIOUR_API_VERSION 1

/*!
 * Macro extracted from :
 * "Why is the new C++ visibility support so useful?"
 * from http://gcc.gnu.org/wiki/Visibility
 */
#if defined _WIN32 || defined _WIN64 || defined __CYGWIN__
#define MGIS_VISIBILITY_IMPORT __declspec(dllimport)
#define MGIS_VISIBILITY_EXPORT __declspec(dllexport)
#define MGIS_VISIBILITY_LOCAL
#else /* defined _WIN32 || defined __CYGWIN__ */
#if (defined __GNUC__) && (!defined __INTEL_COMPILER) && \
    (!defined __NVCOMPILER) && (!defined __clang__)
#if __GNUC__ >= 4
#define MGIS_VISIBILITY_IMPORT [[gnu::visibility("default")]]
#define MGIS_VISIBILITY_EXPORT [[gnu::visibility("default")]]
#define MGIS_VISIBILITY_LOCAL [[gnu::visibility("hidden")]]
#else /*__GNUC__ >= 4 */
#define MGIS_VISIBILITY_IMPORT
#define MGIS_VISIBILITY_EXPORT
#define MGIS_VISIBILITY_LOCAL
#endif /* LIB_MGIS_CONFIG_HXX */
#elif defined __INTEL_COMPILER
#define MGIS_VISIBILITY_IMPORT [[gnu::visibility("default")]]
#define MGIS_VISIBILITY_EXPORT [[gnu::visibility("default")]]
#define MGIS_VISIBILITY_LOCAL [[gnu::visibility("hidden")]]
#elif (defined __NVCOMPILER)
#define MGIS_VISIBILITY_IMPORT __attribute__((visibility("default")))
#define MGIS_VISIBILITY_EXPORT __attribute__((visibility("default")))
#define MGIS_VISIBILITY_LOCAL __attribute__((visibility("hidden")))
#elif defined __clang__
#if __clang_major__ >= 18
#define MGIS_VISIBILITY_IMPORT [[gnu::visibility("default")]]
#define MGIS_VISIBILITY_EXPORT [[gnu::visibility("default")]]
#define MGIS_VISIBILITY_LOCAL [[gnu::visibility("hidden")]]
#else /* __clang_major__ >= 18 */
#define MGIS_VISIBILITY_IMPORT __attribute__((visibility("default")))
#define MGIS_VISIBILITY_EXPORT __attribute__((visibility("default")))
#define MGIS_VISIBILITY_LOCAL __attribute__((visibility("hidden")))
#endif /* __clang_major__ >= 18 */
#else
#define MGIS_VISIBILITY_IMPORT
#define MGIS_VISIBILITY_EXPORT
#define MGIS_VISIBILITY_LOCAL
#endif /* LIB_MGIS_CONFIG_HXX */
#endif /* LIB_MGIS_CONFIG_HXX */

#ifdef MGIS_REAL_TYPE
/*! \brief alias to the numeric type used in the library */
typedef MGIS_REAL_TYPE mgis_real;
#else  /* MGIS_REAL_TYPE */
/*! \brief alias to the numeric type used in the library */
typedef double mgis_real;
#endif /* MGIS_REAL_TYPE */

//! alias to the index type type used in the library
typedef size_t mgis_size_type;

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

/*!
 * \def MGIS_PP_JOIN
 * \brief this macro joins joins its two arguments together.
 *
 * This macro was taken from the boost library:
 * - http://boost.org/
 *
 * The following piece of macro magic joins the two
 * arguments together, even when one of the arguments is
 * itself a macro (see 16.3.1 in C++ standard).  The key
 * is that macro expansion of macro arguments does not
 * occur in MGIS_PP_DO_JOIN2 but does in MGIS_PP_DO_JOIN.
 */
#define MGIS_PP_JOIN(X, Y) MGIS_PP_DO_JOIN(X, Y)

/*!
 * \def MGIS_PP_DO_JOIN
 * \brief An helper macro for MGIS_PP_JOIN
 * \see MGIS_PP_JOIN
 */
#define MGIS_PP_DO_JOIN(X, Y) MGIS_PP_DO_JOIN2(X, Y)

/*!
 * \def MGIS_PP_DO_JOIN2
 * \brief An helper macro for MGIS_PP_JOIN
 * \see MGIS_PP_JOIN
 */
#define MGIS_PP_DO_JOIN2(X, Y) X##Y

/*!
 * \brief macro defining a unique local variable name for a temporary variable
 * that it not meant to be seen by the end-user
 */
#define MGIS_TEMPORARY_VARIABLE(X) \
  MGIS_PP_JOIN(MGIS_PP_JOIN(MGIS_PP_JOIN(mgis_temporary_, X), _), __LINE__)

#endif /* LIB_MGIS_CONFIG_C_H */
