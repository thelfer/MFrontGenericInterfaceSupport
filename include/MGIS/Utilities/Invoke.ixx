/*!
 * \file   MGIS/Utilities/Invoke.ixx
 * \brief  This header implements the invoke function to call function objects
 * that may throw exceptions \date   07/11/2022
 */

#ifndef LIB_MGIS_UTILITIES_INVOKE_IXX
#define LIB_MGIS_UTILITIES_INVOKE_IXX

#include <cerrno>
#include <utility>
#include <cstring>

namespace mgis {

#ifdef MANTA_USE_SOURCE_LOCATION_INFORMATION

  template <typename F, typename... ArgumentsTypes>
  invoke_result_t<F, ArgumentsTypes...> invokeCheckErrno(
      Context &ctx,
      const std::source_location &l,
      F &&f,
      ArgumentsTypes &&...args) noexcept {
    constexpr auto is_void =
        std::is_void_v<std::invoke_result_t<F, ArgumentsTypes...>>;
    if constexpr (is_void) {
      try {
        const auto errno_old = errno;
        errno = 0;
        std::invoke(std::forward<F>(f), std::forward<ArgumentsTypes>(args)...);
        const auto local_errno = errno;
        errno = errno_old;
        if (local_errno != 0) {
          return ctx.registerErrorMessage(std::string{strerror(local_errno)},
                                          l);
        }
      } catch (...) {
        return registerExceptionInErrorBacktrace(ctx, l);
      }
      return true;
    } else {
      try {
        const auto errno_old = errno;
        errno = 0;
        auto r = invoke_result_t<F, ArgumentsTypes...>{};
        r = std::invoke(std::forward<F>(f),
                        std::forward<ArgumentsTypes>(args)...);
        const auto local_errno = errno;
        errno = errno_old;
        if (local_errno != 0) {
          return ctx.registerErrorMessage(std::string{strerror(local_errno)},
                                          l);
        }
        return r;
      } catch (...) {
        registerExceptionInErrorBacktrace(ctx, l);
      }
      return {};
    }
  }

  template <typename F, typename... ArgumentsTypes>
  invoke_result_t<F, ArgumentsTypes...> invoke(
      Context &ctx,
      const std::source_location &l,
      F &&f,
      ArgumentsTypes &&...args) noexcept {
    constexpr auto is_void =
        std::is_void_v<std::invoke_result_t<F, ArgumentsTypes...>>;
    if constexpr (is_void) {
      try {
        std::invoke(std::forward<F>(f), std::forward<ArgumentsTypes>(args)...);
      } catch (...) {
        return registerExceptionInErrorBacktrace(ctx, l);
      }
      return true;
    } else {
      try {
        auto r = invoke_result_t<F, ArgumentsTypes...>{};
        r = std::invoke(std::forward<F>(f),
                        std::forward<ArgumentsTypes>(args)...);
        return r;
      } catch (...) {
        registerExceptionInErrorBacktrace(ctx, l);
      }
      return {};
    }
  }  // end of invoke

#endif

  template <typename F, typename... ArgumentsTypes>
  invoke_result_t<F, ArgumentsTypes...> invokeCheckErrno(
      Context &ctx, F &&f, ArgumentsTypes &&...args) noexcept {
    constexpr auto is_void =
        std::is_void_v<std::invoke_result_t<F, ArgumentsTypes...>>;
    if constexpr (is_void) {
      try {
        const auto errno_old = errno;
        errno = 0;
        std::invoke(std::forward<F>(f), std::forward<ArgumentsTypes>(args)...);
        const auto local_errno = errno;
        errno = errno_old;
        if (local_errno != 0) {
          return ctx.registerErrorMessageWithoutSourceLocation(
              std::string{strerror(local_errno)});
        }
      } catch (...) {
        return registerExceptionInErrorBacktraceWithoutSourceLocation(ctx);
      }
      return true;
    } else {
      try {
        const auto errno_old = errno;
        errno = 0;
        auto r = invoke_result_t<F, ArgumentsTypes...>{};
        r = std::invoke(std::forward<F>(f),
                        std::forward<ArgumentsTypes>(args)...);
        const auto local_errno = errno;
        errno = errno_old;
        if (local_errno != 0) {
          return ctx.registerErrorMessageWithoutSourceLocation(
              std::string{strerror(local_errno)});
        }
        return r;
      } catch (...) {
        std::ignore =
            registerExceptionInErrorBacktraceWithoutSourceLocation(ctx);
      }
      return {};
    }
  }

  template <typename F, typename... ArgumentsTypes>
  invoke_result_t<F, ArgumentsTypes...> invoke(
      Context &ctx, F &&f, ArgumentsTypes &&...args) noexcept {
    constexpr auto is_void =
        std::is_void_v<std::invoke_result_t<F, ArgumentsTypes...>>;
    if constexpr (is_void) {
      try {
        std::invoke(std::forward<F>(f), std::forward<ArgumentsTypes>(args)...);
      } catch (...) {
        return registerExceptionInErrorBacktraceWithoutSourceLocation(ctx);
      }
      return true;
    } else {
      try {
        auto r = invoke_result_t<F, ArgumentsTypes...>{};
        r = std::invoke(std::forward<F>(f),
                        std::forward<ArgumentsTypes>(args)...);
        return r;
      } catch (...) {
        std::ignore =
            registerExceptionInErrorBacktraceWithoutSourceLocation(ctx);
      }
      return {};
    }
  }  // end of invoke

}  // end of namespace mgis

#endif /* LIB_MGIS_UTILITIES_INVOKE_IXX */