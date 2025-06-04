/*!
 * \file   MGIS/AbstractErrorHandler.hxx
 * \brief
 * \author Thomas Helfer
 * \date   04/06/2025
 */

#ifndef LIB_MGIS_ABSTRACTERRORHANDLER_HXX
#define LIB_MGIS_ABSTRACTERRORHANDLER_HXX

#include "MGIS/Config.hxx"
#include "MGIS/InvalidResult.hxx"

// source_location  seems badly supported by the current compilers (2022):
// - with some compilers, this header is missing
// - with clang this header exists, but std::source_location is undefined
// The following tests tries to handle all those corner cases. In the
// future, this shall be removed...
#ifdef MGIS_DEBUG
#if __has_include(<source_location>)
#ifndef __clang__
#include <source_location>
#include <cstdint>
#define MGIS_USE_SOURCE_LOCATION_INFORMATION
#endif /* __clang */
#else
#endif /* __has_include(<source_location>) */
#endif /* MGIS_DEBUG */

namespace mgis {

  /*!
   *
   */
  struct MGIS_EXPORT AbstractErrorHandler {
#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
    /*!
     * \brief register a new error message
     * \param[in] e: error message
     * \param[in] l: description of the call site
     * \note for convenience, this method always return `false`
     */
    virtual InvalidResult registerErrorMessage(
        const char *const,
        const std::source_location & =
            std::source_location::current()) noexcept;
#else
    /*!
     * \brief register a new error message
     * \param[in] e: error message
     * \note for convenience, this method always return `InvalidResult`
     */
    virtual InvalidResult registerErrorMessage(const char *const) = 0;
#endif
    virtual constexpr ~AbstractErrorHandler() = default;
  };

}  // end of namespace mgis

#endif /* LIB_MGIS_ABSTRACTERRORHANDLER_HXX */
