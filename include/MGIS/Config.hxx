/*!
 * \file   include/MGIS/Config.hxx
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

#ifndef LIB_MGIS_CONFIG_HXX
#define LIB_MGIS_CONFIG_HXX

#include <limits>
#include "MGIS/Config-c.h"

namespace mgis::config {

  /*!
   * \brief policy used in case of contract violation
   */
  enum struct ContractViolationPolicy { ABORT, RAISE };

  /*!
   * \brief boolean variable stating which policy is usd for
   * reporting a contract violation
   *
   * By default, `std::abort` is called.
   *
   * If `MGIS_USE_EXCEPTIONS_FOR_CONTRACT_VIOLATION` is defined
   * `mgis::raise` is called.
   */
  inline constexpr auto contract_violation_policy =
#ifdef MGIS_USE_EXCEPTIONS_FOR_CONTRACT_VIOLATION
      ContractViolationPolicy::RAISE;
#else
      ContractViolationPolicy::ABORT;
#endif

  /*!
   * \brief  error reporting policy used by default when a Context is
   * available.
   */
  enum struct ErrorReportPolicy { INVALIDRESULT, RAISE, ABORT };

  /*!
   * \brief boolean variable stating which policy is usd for
   * reporting a runtime error
   *
   * By default, error are reported by returning an invalid result and
   * registring an error message in a context
   *
   * If MGIS_USE_EXCEPTIONS_FOR_ERROR_REPORTING is defined, registring
   * an error message in a context, `mgis::raise` is called
   *
   * If MGIS_USE_ABORT_FOR_ERROR_REPORTING is defined, registring
   * an error message in a context, the error message is printed
   * an the standard error stream and `std::abort` is called.
   */
  inline constexpr auto error_report_policy =
#ifdef MGIS_USE_EXCEPTIONS_FOR_ERROR_REPORTING
#ifdef MGIS_USE_ABORT_FOR_ERROR_REPORTING
#error \
    "MGIS_USE_EXCEPTIONS_FOR_ERROR_REPORTING and " \
        "MGIS_USE_ABORT_FOR_ERROR_REPORTING " \
        "are mutually exclusive"
#endif /* MGIS_USE_ABORT_FOR_ERROR_REPORTING */
      ErrorReportPolicy::RAISE;
#else /* MGIS_USE_EXCEPTIONS_FOR_ERROR_REPORTING */
#ifdef MGIS_USE_ABORT_FOR_ERROR_REPORTING
      ErrorReportPolicy::ABORT;
#else
      ErrorReportPolicy::INVALIDRESULT;
#endif
#endif /* MGIS_USE_EXCEPTIONS_FOR_ERROR_REPORTING */

}  // end of namespace mgis::config

namespace mgis::attributes {

  /*!
   * \brief an attribute use to indicate that a method or a function as being
   * unsafe without precautions
   */
  struct UnsafeAttribute {};
  /*!
   * \brief an attribute use to indicate that a method may throw or not
   */
  template <bool>
  struct ThrowingAttribute {};

}  // namespace mgis::attributes

namespace mgis {

  //
  inline constexpr auto unsafe = attributes::UnsafeAttribute{};
  //
  inline constexpr auto throwing = attributes::ThrowingAttribute<true>{};
  //
  inline constexpr auto not_throwing = attributes::ThrowingAttribute<false>{};

  //! \brief a simple alias to the the default indexing type used by mgis
  using size_type = mgis_size_type;

  //! \brief alias to the numeric type used
  using real = mgis_real;

  //! \brief a constant whose role is similar to std::dynamic_extent
  inline constexpr size_type dynamic_extent =
      std::numeric_limits<size_type>::max();

}  // end of namespace mgis

#endif /* LIB_MGIS_CONFIG_HXX */
