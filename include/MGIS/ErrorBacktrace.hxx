/*!
 * \file   MGIS/ErrorBacktrace.hxx
 * \brief  This file declares the `ErrorBacktrace` class.
 * \date   04/11/2022
 */

#ifndef LIB_MGIS_ERRORBACKTRACE_HXX
#define LIB_MGIS_ERRORBACKTRACE_HXX 1

#include <string>
#include <memory>
#include <vector>
#include <variant>
#include <utility>
#include <type_traits>
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
   * \brief a structure to keep trace of error messages.
   *
   * The design of this class is inspired by the many solutions offered
   * by the standard library (exceptions, `std::error_code`, etc...)
   */
  class MGIS_EXPORT ErrorBacktrace {
   public:
    //! brief a simple alias
    using ErrorReportFunction = std::string (*)(const int);
    /*!
     * \brief a simple alias listing the main way to report an error.
     *
     * The alternatives way of reporting an error are:
     *
     * - a plain message in the form of a C-string
     * - a message in the form of a std::string
     * - an error code associated with a function allowing to retrieve an
     *   appropriate error message.
     *
     * Building complex error message shall be avoided as much as possible in
     * template code.
     */
    using ErrorReport = std::
        variant<const char *, std::pair<int, ErrorReportFunction>, std::string>;
    //! \brief specify if error reporting shall be fatal
    static void setErrorReportingAsFatal() noexcept;
    //! \brief specify if error reporting shall be fatal
    static void unsetErrorReportingAsFatal() noexcept;
#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
    /*!
     * \brief register a new error message
     * \param[in] e: error code
     * \param[in] l: description of the call site
     * \note for convenience, this method always return `false`
     */
    InvalidResult registerErrorMessage(
        const ErrorReport,
        const std::source_location & =
            std::source_location::current()) noexcept;
#else
    /*!
     * \brief register a new error message
     * \param[in] e: error code
     * \note for convenience, this method always return `InvalidResult`
     */
    InvalidResult registerErrorMessage(const ErrorReport) noexcept;
#endif
    /*!
     * \brief register a new error message
     * \param[in] e: error code
     * \note for convenience, this method always return `InvalidResult`
     */
    InvalidResult registerErrorMessageWithoutSourceLocation(
        const ErrorReport) noexcept;
    //! \brief remove the error messages
    void clearErrorMessages() noexcept;
    /*!
     * \return the registered error messages
     * \param[in] b: boolean. If true, the error message are cleared
     */
    std::string getErrorMessage(const bool = true) noexcept;
    /*!
     * \return the error message without the decoration (file name, function
     * name, line number) \param[in] b: boolean. If true, the error message are
     * cleared
     */
    std::string getRawErrorMessage(const bool = true) noexcept;
    //! \return if the list registered error messages is empty
    bool empty() const noexcept;
    //! \brief destructor
    ~ErrorBacktrace() noexcept;

   private:
    //! structure describing an error message
    struct ErrorMessage {
#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
      //! \brief file name
      const char *file_name = nullptr;
      //! \brief function name
      const char *function_name = nullptr;
      //! \brief line number
      std::uint_least32_t line = 0u;
#endif
      //! \brief structure used to retrieve the error message
      ErrorReport e;
    };
    /*!
     * \return the error message
     * \param[in] e: error report
     */
    static std::string getErrorMessage(const ErrorReport &) noexcept;
    //! \return the registered error messages
    std::string getErrorMessage_() const noexcept;
    //! \return the error message without the decoration (file name, function
    //! name, line number)
    std::string getRawErrorMessage_() const noexcept;
    //! \brief treat the case when error reporting is fatal
    void treatFatalCase_() const noexcept;
    //! \brief list of registered error message
    std::vector<ErrorMessage> errorMessages_;
  };  // end of ErrorBacktrace

#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
  /*!
   * \brief a custom Lippincott-like function that extract error messages from
   * an exception an \param[out] e: error back trace handler \param[in] l:
   * description of the call site
   */
  MGIS_EXPORT InvalidResult registerExceptionInErrorBacktrace(
      ErrorBacktrace &,
      const std::source_location & = std::source_location::current()) noexcept;
#else
  /*!
   * \brief a custom Lippincott-like function that extract error messages from
   * an exception an \param[out] e: error back trace handler
   */
  MGIS_EXPORT InvalidResult
  registerExceptionInErrorBacktrace(ErrorBacktrace &) noexcept;
#endif
  /*!
   * \brief a custom Lippincott-like function that extract error messages from
   * an exception an \param[out] e: error back trace handler
   */
  MGIS_EXPORT InvalidResult
  registerExceptionInErrorBacktraceWithoutSourceLocation(
      ErrorBacktrace &) noexcept;

}  // end of namespace mgis

#endif /* LIB_MGIS_ERRORBACKTRACE_HXX */