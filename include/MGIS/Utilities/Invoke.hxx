/*!
 * \file   MGIS/Utilities/Invoke.hxx
 * \brief  This header declares the invoke function to call function objects
 * that may throw exceptions \date   07/11/2022
 */

#ifndef LIB_MGIS_UTILITIES_INVOKE_HXX
#define LIB_MGIS_UTILITIES_INVOKE_HXX

#include <optional>
#include <functional>
#include <type_traits>
#include "MGIS/Config.hxx"
#include "MGIS/InvalidResult.hxx"
#include "MGIS/Context.hxx"

/*!
 * \def MGIS_INVOKE(ctx, ...)
 * \brief a simple wrapper around the `invoke` function which handles the
 * declaration of the source location if needed
 *
 * \param[in, out] ctx: execution context
 *
 * \warning If the function invoked has optional arguments (with default
 * values), all the optional arguments must be provided, else the code will not
 * compile with an obscure error message \note a version check the errno value
 * while the other don't (names are self descriptive)
 */
#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
#define MGIS_INVOKE(e, ...) \
  ::mgis::invoke(e, std::source_location::current(), __VA_ARGS__)
#define MGIS_INVOKE_CHECK_ERRNO(e, ...) \
  ::mgis::invokeCheckErrno(e, std::source_location::current(), __VA_ARGS__)
#else
#define MGIS_INVOKE(e, ...) ::mgis::invoke(e, __VA_ARGS__)
#define MGIS_INVOKE_CHECK_ERRNO(e, ...) ::mgis::invokeCheckErrno(e, __VA_ARGS__)
#endif

/*!
 * \def MGIS_TRY_INVOKE_VOID(ctx, ...)
 * \brief a simple wrapper around the `invoke` function when the returned type
 * of the callable is void
 *
 * \param[in, out] ctx: execution context
 *
 * \warning If the
 * function invoked has optional arguments (with default values), all the
 * optional arguments must be provided, else the code will not compile with an
 * obscure error message \note a version check the errno value while the other
 * don't (names are self descriptive)
 */
#define MGIS_TRY_INVOKE_VOID(e, ...)  \
  if (!MGIS_INVOKE(e, __VA_ARGS__)) { \
    return InvalidResult{};           \
  }
#define MGIS_TRY_INVOKE_VOID_CHECK_ERRNO(e, ...)  \
  if (!MGIS_INVOKE_CHECK_ERRNO(e, __VA_ARGS__)) { \
    return InvalidResult{};                       \
  }

/*!
 * \def MGIS_TRY_INVOKE(Type, e, ...)
 * \brief a simple wrapper around the `invoke` function
 *
 * \param[inout] ctx: execution context
 * \param[out] v: variable name which will contain the returned value of the
 * function invoked. This variable is declared inside the macro and so must not
 * be declared before.
 *
 * \warning If the function invoked has optional arguments (with default
 * values), all the optional arguments must be provided, else the code will not
 * compile with an obscure error message
 * \note a version check the errno value while the other don't (names are self
 * descriptive)
 */
#define MGIS_TRY_INVOKE(v, e, ...)                                     \
  const auto MGIS_TEMPORARY_VARIABLE(v) = MGIS_INVOKE(e, __VA_ARGS__); \
  if (!MGIS_TEMPORARY_VARIABLE(v).hasValue()) {                        \
    return InvalidResult{};                                            \
  }                                                                    \
  auto &v = *(MGIS_TEMPORARY_VARIABLE(v));
#define MGIS_TRY_INVOKE_CHECK_ERRNO(v, e, ...)  \
  const auto MGIS_TEMPORARY_VARIABLE(v) =       \
      MGIS_INVOKE_CHECK_ERRNO(e, __VA_ARGS__);  \
  if (!MGIS_TEMPORARY_VARIABLE(v).hasValue()) { \
    return InvalidResult{};                     \
  }                                             \
  auto &v = *(MGIS_TEMPORARY_VARIABLE(v));

namespace mgis {

  //! \brief a simple alias
  template <typename F, typename... ArgumentsTypes>
  using invoke_result_t = std::conditional_t<
      std::is_void_v<std::invoke_result_t<F, ArgumentsTypes...>>,
      bool,
      std::optional<std::invoke_result_t<F, ArgumentsTypes...>>>;

#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION

  /*!
   * \brief invoke the given callable object
   * param[inout] ctx: execution context
   * \param[in] l: source location
   * \param[in] F: callable object
   * \param[in] args: arguments
   * \warning If the function invoked has optional arguments (with default
   * values), all the optional arguments must be provided, else the code will
   * not compile with an obscure error message \note a version check the errno
   * value while the other don't (names are self descriptive)
   */
  template <typename F, typename... ArgumentsTypes>
  invoke_result_t<F, ArgumentsTypes...> invoke(Context &,
                                               const std::source_location &,
                                               F &&,
                                               ArgumentsTypes &&...) noexcept;
  template <typename F, typename... ArgumentsTypes>
  invoke_result_t<F, ArgumentsTypes...> invokeCheckErrno(
      Context &,
      const std::source_location &,
      F &&,
      ArgumentsTypes &&...) noexcept;

#endif /* MGIS_USE_SOURCE_LOCATION_INFORMATION */

  /*!
   * \brief invoke the given callable object
   * param[inout] ctx: execution context
   * \param[in] F: callable object
   * \param[in] args: arguments
   * \warning If the function invoked has optional arguments (with default
   * values), all the optional arguments must be provided, else the code will
   * not compile with an obscure error message \note a version check the errno
   * value while the other don't (names are self descriptive)
   */
  template <typename F, typename... ArgumentsTypes>
  invoke_result_t<F, ArgumentsTypes...> invoke(Context &,
                                               F &&,
                                               ArgumentsTypes &&...) noexcept;
  template <typename F, typename... ArgumentsTypes>
  invoke_result_t<F, ArgumentsTypes...> invokeCheckErrno(
      Context &, F &&, ArgumentsTypes &&...) noexcept;

}  // end of namespace mgis


#include "MGIS/Utilities/Invoke.ixx"

#endif /* LIB_MGIS_UTILITIES_INVOKE_HXX */
