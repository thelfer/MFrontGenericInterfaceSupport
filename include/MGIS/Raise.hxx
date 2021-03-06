/*!
 * \file   include/MGIS/Raise.hxx
 * \brief  declaration of the `raise` function.
 * \author Thomas Helfer
 * \date   20/06/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_RAISE_HXX
#define LIB_MGIS_RAISE_HXX

#include <utility>
#include <stdexcept>
#include "MGIS/Config.hxx"

namespace mgis {

  //! \brief a simple alias
  using ExceptionHandler = void (*)(void);

  /*!
   * \brief set an exception handler
   * \param[in] h: exception handler
   *
   * \code{.cpp}
   * void handler() {
   *   try {
   *     throw;
   *   } catch (std::exception& e) {
   *     std::cerr << e.what() << '\n';
   *   } catch (...) {
   *     std::cerr << "unknown exception thrown";
   *   }
   *   std::abort();
   * }  // end of handler
   *
   * mgis::setExceptionHandler(handler);
   * mgis::raise<std::logic_error>("something went wrong");
   * \endcode
   */
  MGIS_EXPORT void setExceptionHandler(ExceptionHandler);

  /*!
   * \brief return a registred exception handler, nullptr if none were
   * registred.
   */
  MGIS_EXPORT ExceptionHandler getExceptionHandler();

  /*!
   * \brief a small wrapper used to build the exception outside the
   * `throw` statement. As most exception's classes constructors may
   * throw, this avoids undefined behaviour as reported by the
   * `cert-err60-cpp` warning of `clang-tidy` (thrown exception type
   * is not nothrow copy constructible).
   * \tparam Exception: type of the exception to be thrown.
   */
  template <typename Exception = std::runtime_error>
  [[noreturn]] MGIS_VISIBILITY_LOCAL inline void raise();

  /*!
   * \brief a small wrapper used to build the exception outside the
   * `throw` statement. As most exception's classes constructors may
   * throw, this avoids undefined behaviour as reported by the
   * `cert-err60-cpp` warning of `clang-tidy` (thrown exception type
   * is not nothrow copy constructible).
   * \tparam Exception: type of the exception to be thrown.
   * \tparam Args: type of the arguments passed to the exception'
   * constructor.
   * \param[in] a: arguments passed to the exception' constructor.
   */
  template <typename Exception = std::runtime_error, typename... Args>
  [[noreturn]] MGIS_VISIBILITY_LOCAL inline void raise(Args&&...);

  /*!
   * \brief raise an exception if the first argument is `true`.
   * \tparam Exception: type of the exception to be thrown.
   * \param[in] b: condition to be checked. If `true`, an exception is
   * thrown.
   */
  template <typename Exception = std::runtime_error>
  MGIS_VISIBILITY_LOCAL inline void raise_if(const bool);

  /*!
   * \brief raise an exception if the first argument is `true`.
   * \tparam Exception: type of the exception to be thrown.
   * \tparam Args: type of the arguments passed to the exception'
   * constructor.
   * \param[in] b: condition to be checked. If `true`, an exception is
   * thrown.
   * \param[in] a: arguments passed to the exception' constructor.
   */
  template <typename Exception = std::runtime_error, typename... Args>
  MGIS_VISIBILITY_LOCAL inline void raise_if(const bool, Args&&...);

}  // end of namespace mgis

#include "MGIS/Raise.ixx"

#endif /* LIB_MGIS_RAISE_HXX */
