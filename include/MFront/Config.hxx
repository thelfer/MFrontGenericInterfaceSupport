/*!
 * \file   include/MFront/Config.hxx
 * \brief    
 * \author th202608
 * \date   19/06/2018
 */

#ifndef LIB_MFRONT_CONFIG_HXX
#define LIB_MFRONT_CONFIG_HXX

#include <cstddef>

/*!
 * Macro extracted from :
 * "Why is the new C++ visibility support so useful?"
 * from http://gcc.gnu.org/wiki/Visibility
 */
#if defined _WIN32 || defined _WIN64 ||defined __CYGWIN__
#define MFRONT_VISIBILITY_IMPORT __declspec(dllimport)
#define MFRONT_VISIBILITY_EXPORT __declspec(dllexport)
#define MFRONT_VISIBILITY_LOCAL
#else /* defined _WIN32 || defined __CYGWIN__ */
#if (defined __GNUC__) && (! defined __INTEL_COMPILER)
#if __GNUC__ >= 4
#define MFRONT_VISIBILITY_IMPORT __attribute__((visibility("default")))
#define MFRONT_VISIBILITY_EXPORT __attribute__((visibility("default")))
#define MFRONT_VISIBILITY_LOCAL  __attribute__((visibility("hidden")))
#else /* __GNUC__ >= 4 */
#define MFRONT_VISIBILITY_IMPORT
#define MFRONT_VISIBILITY_EXPORT
#define MFRONT_VISIBILITY_LOCAL
#endif /* __GNUC__ >= 4 */
#elif defined __INTEL_COMPILER
#define MFRONT_VISIBILITY_IMPORT __attribute__((visibility("default")))
#define MFRONT_VISIBILITY_EXPORT __attribute__((visibility("default")))
#define MFRONT_VISIBILITY_LOCAL  __attribute__((visibility("hidden")))
#else /* defined __INTEL_COMPILER */
#define MFRONT_VISIBILITY_IMPORT
#define MFRONT_VISIBILITY_EXPORT
#define MFRONT_VISIBILITY_LOCAL
#endif /* defined __INTEL_COMPILER */
#endif /* defined _WIN32 || defined _WIN64 ||defined __CYGWIN__ */

#if defined _WIN32 || defined _WIN64 || defined __CYGWIN__
#if defined MFrontGenericInterface_EXPORTS
#define MFRONT_EXPORT MFRONT_VISIBILITY_EXPORT
#else
#ifndef MFRONT_STATIC_BUILD
#define MFRONT_EXPORT MFRONT_VISIBILITY_IMPORT
#else
#define MFRONT_EXPORT
#endif
#endif
#else
#define MFRONT_EXPORT MFRONT_VISIBILITY_EXPORT
#endif /* LIB_MFRONT_CONFIG_HXX */

#define MFRONT_INLINE inline
#define MFRONT_NORETURN [[noreturn]]

namespace mfront {

  namespace behaviour {
  
    //! a simple alias
    using size_type = size_t;

#ifdef MFRONT_REAL_TYPE
    //! alias to the numeric type used
    using real = MFRONT_REAL_TYPE;
#else
    //! alias to the numeric type used
    using real = double;
#endif

  }  // end of namespace behaviour

}  // end of namespace mfront

#endif /* LIB_MFRONT_CONFIG_HXX */
