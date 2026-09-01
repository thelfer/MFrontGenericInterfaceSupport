/*!
 * \file   MIGS/Utilities/Construct.ixx
 * \brief  This file implements the template functions declared
 *         in the `MGIS/Construct.hx` header.
 * \date   04/11/2022
 */

#ifndef LIB_MGIS_UTILITIES_CONSTRUCT_IXX
#define LIB_MGIS_UTILITIES_CONSTRUCT_IXX

namespace mgis::internals {

#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION

  template <typename Type, typename... ArgumentsTypes>
  std::optional<Type> construct_impl(Context &ctx,
                                     const std::source_location &l,
                                     ArgumentsTypes &&...args) noexcept requires
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...> {
    using NonConstType = std::remove_const_t<Type>;
    if constexpr (std::is_nothrow_constructible_v<NonConstType,
                                                  ArgumentsTypes...>) {
      return std::make_optional<Type>(std::forward<ArgumentsTypes>(args)...);
    } else {
      try {
        return std::make_optional<Type>(std::forward<ArgumentsTypes>(args)...);
      } catch (...) {
        registerExceptionInErrorBacktrace(ctx, l);
      }
      return {};
    }
  }  // end of construct

  template <typename Type, typename... ArgumentsTypes>
  std::unique_ptr<Type> make_unique_impl(
      Context &ctx,
      const std::source_location &l,
      ArgumentsTypes &&...args) noexcept requires
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...> {
    using NonConstType = std::remove_const_t<Type>;
    if constexpr (std::is_nothrow_constructible_v<NonConstType,
                                                  ArgumentsTypes...>) {
      return std::make_unique<Type>(std::forward<ArgumentsTypes>(args)...);
    } else {
      try {
        return std::make_unique<Type>(std::forward<ArgumentsTypes>(args)...);
      } catch (...) {
        registerExceptionInErrorBacktrace(ctx, l);
      }
      return {};
    }
  }  // end of make_unique

  template <typename BaseType, typename Type, typename... ArgumentsTypes>
  std::unique_ptr<BaseType> make_unique_as_impl(
      Context &ctx,
      const std::source_location &l,
      ArgumentsTypes &&...args) noexcept requires
      std::is_base_of_v<BaseType, Type> &&
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...> {
    using NonConstType = std::remove_const_t<Type>;
    if constexpr (std::is_nothrow_constructible_v<NonConstType,
                                                  ArgumentsTypes...>) {
      auto *const p = new NonConstType{std::forward<ArgumentsTypes>(args)...};
      return std::unique_ptr<BaseType>{p};
    } else {
      try {
        auto *const p = new NonConstType{std::forward<ArgumentsTypes>(args)...};
        return std::unique_ptr<BaseType>{p};
      } catch (...) {
        registerExceptionInErrorBacktrace(ctx, l);
      }
      return {};
    }
  }  // end of make_unique_as

  template <typename Type, typename... ArgumentsTypes>
  std::shared_ptr<Type> make_shared_impl(
      Context &ctx,
      const std::source_location &l,
      ArgumentsTypes &&...args) noexcept requires
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...> {
    using NonConstType = std::remove_const_t<Type>;
    if constexpr (std::is_nothrow_constructible_v<NonConstType,
                                                  ArgumentsTypes...>) {
      return std::make_shared<Type>(std::forward<ArgumentsTypes>(args)...);
    } else {
      try {
        return std::make_shared<Type>(std::forward<ArgumentsTypes>(args)...);
      } catch (...) {
        registerExceptionInErrorBacktrace(ctx, l);
      }
      return {};
    }
  }  // end of make_shared

  template <typename BaseType, typename Type, typename... ArgumentsTypes>
  std::shared_ptr<BaseType> make_shared_as_impl(
      Context &ctx,
      const std::source_location &l,
      ArgumentsTypes &&...args) noexcept requires
      std::is_base_of_v<BaseType, Type> &&
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...> {
    using NonConstType = std::remove_const_t<Type>;
    if constexpr (std::is_nothrow_constructible_v<NonConstType,
                                                  ArgumentsTypes...>) {
      auto *const p = new NonConstType{std::forward<ArgumentsTypes>(args)...};
      return std::shared_ptr<BaseType>{p};
    } else {
      try {
        auto *const p = new NonConstType{std::forward<ArgumentsTypes>(args)...};
        return std::shared_ptr<BaseType>{p};
      } catch (...) {
        registerExceptionInErrorBacktrace(ctx, l);
      }
      return {};
    }
  }  // end of make_shared_as

#endif /* MGIS_USE_SOURCE_LOCATION_INFORMATION */

  template <typename Type, typename... ArgumentsTypes>
  std::optional<Type> construct_impl(Context &ctx,
                                     ArgumentsTypes &&...args) noexcept requires
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...> {
    using NonConstType = std::remove_const_t<Type>;
    if constexpr (std::is_nothrow_constructible_v<NonConstType,
                                                  ArgumentsTypes...>) {
      return Type{std::forward<ArgumentsTypes>(args)...};
    } else {
      try {
        return std::make_optional<Type>(std::forward<ArgumentsTypes>(args)...);
      } catch (...) {
        std::ignore =
            registerExceptionInErrorBacktraceWithoutSourceLocation(ctx);
      }
      return {};
    }
  }  // end of construct

  template <typename Type, typename... ArgumentsTypes>
  std::unique_ptr<Type> make_unique_impl(
      Context &ctx, ArgumentsTypes &&...args) noexcept requires
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...> {
    using NonConstType = std::remove_const_t<Type>;
    if constexpr (std::is_nothrow_constructible_v<NonConstType,
                                                  ArgumentsTypes...>) {
      return std::make_unique<Type>(std::forward<ArgumentsTypes>(args)...);
    } else {
      try {
        return std::make_unique<Type>(std::forward<ArgumentsTypes>(args)...);
      } catch (...) {
        std::ignore =
            registerExceptionInErrorBacktraceWithoutSourceLocation(ctx);
      }
      return {};
    }
  }  // end of make_unique

  template <typename BaseType, typename Type, typename... ArgumentsTypes>
  std::unique_ptr<BaseType> make_unique_as_impl(
      Context &ctx, ArgumentsTypes &&...args) noexcept requires
      std::is_base_of_v<BaseType, Type> &&
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...> {
    using NonConstType = std::remove_const_t<Type>;
    if constexpr (std::is_nothrow_constructible_v<NonConstType,
                                                  ArgumentsTypes...>) {
      auto *const p = new NonConstType{std::forward<ArgumentsTypes>(args)...};
      return std::unique_ptr<BaseType>{p};
    } else {
      try {
        auto *const p = new NonConstType{std::forward<ArgumentsTypes>(args)...};
        return std::unique_ptr<BaseType>{p};
      } catch (...) {
        std::ignore =
            registerExceptionInErrorBacktraceWithoutSourceLocation(ctx);
      }
      return {};
    }
  }  // end of make_unique_as

  template <typename Type, typename... ArgumentsTypes>
  std::shared_ptr<Type> make_shared_impl(
      Context &ctx, ArgumentsTypes &&...args) noexcept requires
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...> {
    using NonConstType = std::remove_const_t<Type>;
    if constexpr (std::is_nothrow_constructible_v<NonConstType,
                                                  ArgumentsTypes...>) {
      return std::make_shared<Type>(std::forward<ArgumentsTypes>(args)...);
    } else {
      try {
        return std::make_shared<Type>(std::forward<ArgumentsTypes>(args)...);
      } catch (...) {
        std::ignore =
            registerExceptionInErrorBacktraceWithoutSourceLocation(ctx);
      }
      return {};
    }
  }  // end of make_shared

  template <typename BaseType, typename Type, typename... ArgumentsTypes>
  std::shared_ptr<BaseType> make_shared_as_impl(
      Context &ctx, ArgumentsTypes &&...args) noexcept requires
      std::is_base_of_v<BaseType, Type> &&
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...> {
    using NonConstType = std::remove_const_t<Type>;
    if constexpr (std::is_nothrow_constructible_v<NonConstType,
                                                  ArgumentsTypes...>) {
      auto *const p = new NonConstType{std::forward<ArgumentsTypes>(args)...};
      return std::shared_ptr<BaseType>{p};
    } else {
      try {
        auto *const p = new NonConstType{std::forward<ArgumentsTypes>(args)...};
        return std::shared_ptr<BaseType>{p};
      } catch (...) {
        std::ignore =
            registerExceptionInErrorBacktraceWithoutSourceLocation(ctx);
      }
      return {};
    }
  }  // end of make_shared_as

}  // namespace mgis::internals

namespace mgis {

#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION

  template <typename Type, typename... ArgumentsTypes>
  std::optional<Type> construct(Context &ctx,
                                const std::source_location &l,
                                ArgumentsTypes &&...args) noexcept
      requires ::mgis::internals::is_constructible<Type, ArgumentsTypes...> {
    if constexpr (::mgis::internals::is_constructible_with_context<
                      Type, ArgumentsTypes...>) {
      return ::mgis::internals::construct_impl<Type>(
          ctx, l, ctx, std::forward<ArgumentsTypes>(args)...);
    } else {
      return ::mgis::internals::construct_impl<Type>(
          ctx, l, ctx, std::forward<ArgumentsTypes>(args)...);
    }
  }  // end of construct

  template <typename Type, typename... ArgumentsTypes>
  std::unique_ptr<Type> make_unique(Context &ctx,
                                    const std::source_location &l,
                                    ArgumentsTypes &&...args) noexcept
      requires ::mgis::internals::is_constructible<Type, ArgumentsTypes...> {
    if constexpr (::mgis::internals::is_constructible_with_context<
                      Type, ArgumentsTypes...>) {
      return ::mgis::internals::make_unique_impl<Type>(
          ctx, l, ctx, std::forward<ArgumentsTypes>(args)...);
    } else {
      return ::mgis::internals::make_unique_impl<Type>(
          ctx, l, std::forward<ArgumentsTypes>(args)...);
    }
  }  // end of make_unique

  template <typename BaseType, typename Type, typename... ArgumentsTypes>
  std::unique_ptr<BaseType> make_unique_as(Context &ctx,
                                           const std::source_location &l,
                                           ArgumentsTypes &&...args) noexcept
      requires std::is_base_of_v<BaseType, Type> &&
      ::mgis::internals::is_constructible<Type, ArgumentsTypes...> {
    if constexpr (::mgis::internals::is_constructible_with_context<
                      Type, ArgumentsTypes...>) {
      return ::mgis::internals::make_unique_imp_asl<BaseType, Type>(
          ctx, l, ctx, std::forward<ArgumentsTypes>(args)...);
    } else {
      return ::mgis::internals::make_unique_imp_asl<BaseType, Type>(
          ctx, l, std::forward<ArgumentsTypes>(args)...);
    }
  }  // end of make_unique_as

  template <typename Type, typename... ArgumentsTypes>
  std::shared_ptr<Type> make_shared(Context &ctx,
                                    const std::source_location &l,
                                    ArgumentsTypes &&...args) noexcept
      requires ::mgis::internals::is_constructible<Type, ArgumentsTypes...> {
    if constexpr (::mgis::internals::is_constructible_with_context<
                      Type, ArgumentsTypes...>) {
      return ::mgis::internals::make_shared_impl<Type>(
          ctx, l, ctx, std::forward<ArgumentsTypes>(args)...);
    } else {
      return ::mgis::internals::make_shared_impl<Type>(
          ctx, l, std::forward<ArgumentsTypes>(args)...);
    }
  }  // end of make_shared

  template <typename BaseType, typename Type, typename... ArgumentsTypes>
  std::shared_ptr<BaseType> make_shared_as(Context &ctx,
                                           const std::source_location &l,
                                           ArgumentsTypes &&...args) noexcept
      requires std::is_base_of_v<BaseType, Type> &&
      ::mgis::internals::is_constructible<Type, ArgumentsTypes...> {
    if constexpr (::mgis::internals::is_constructible_with_context<
                      Type, ArgumentsTypes...>) {
      return ::mgis::internals::make_shared_as_impl<BaseType, Type>(
          ctx, l, ctx, std::forward<ArgumentsTypes>(args)...);
    } else {
      return ::mgis::internals::make_shared_as_impl<BaseType, Type>(
          ctx, l, std::forward<ArgumentsTypes>(args)...);
    }
  }  // end of make_shared_as

#endif /* MGIS_USE_SOURCE_LOCATION_INFORMATION */

  template <typename Type, typename... ArgumentsTypes>
  std::optional<Type> construct(Context &ctx, ArgumentsTypes &&...args) noexcept
      requires ::mgis::internals::is_constructible<Type, ArgumentsTypes...> {
    if constexpr (::mgis::internals::is_constructible_with_context<
                      Type, ArgumentsTypes...>) {
      return ::mgis::internals::construct_impl<Type>(
          ctx, ctx, std::forward<ArgumentsTypes>(args)...);
    } else {
      return ::mgis::internals::construct_impl<Type>(
          ctx, std::forward<ArgumentsTypes>(args)...);
    }
  }  // end of construct

  template <typename Type, typename... ArgumentsTypes>
  std::unique_ptr<Type> make_unique(Context &ctx,
                                    ArgumentsTypes &&...args) noexcept
      requires ::mgis::internals::is_constructible<Type, ArgumentsTypes...> {
    if constexpr (::mgis::internals::is_constructible_with_context<
                      Type, ArgumentsTypes...>) {
      return ::mgis::internals::make_unique_impl<Type>(
          ctx, ctx, std::forward<ArgumentsTypes>(args)...);
    } else {
      return ::mgis::internals::make_unique_impl<Type>(
          ctx, std::forward<ArgumentsTypes>(args)...);
    }
  }  // end of make_unique

  template <typename BaseType, typename Type, typename... ArgumentsTypes>
  std::unique_ptr<BaseType> make_unique_as(Context &ctx,
                                           ArgumentsTypes &&...args) noexcept
      requires std::is_base_of_v<BaseType, Type> &&
      ::mgis::internals::is_constructible<Type, ArgumentsTypes...> {
    if constexpr (::mgis::internals::is_constructible_with_context<
                      Type, ArgumentsTypes...>) {
      return ::mgis::internals::make_unique_as_impl<BaseType, Type>(
          ctx, ctx, std::forward<ArgumentsTypes>(args)...);
    } else {
      return ::mgis::internals::make_unique_as_impl<BaseType, Type>(
          ctx, std::forward<ArgumentsTypes>(args)...);
    }
  }  // end of make_unique_as

  template <typename Type, typename... ArgumentsTypes>
  std::shared_ptr<Type> make_shared(Context &ctx,
                                    ArgumentsTypes &&...args) noexcept
      requires ::mgis::internals::is_constructible<Type, ArgumentsTypes...> {
    if constexpr (::mgis::internals::is_constructible_with_context<
                      Type, ArgumentsTypes...>) {
      return ::mgis::internals::make_shared_impl<Type>(
          ctx, ctx, std::forward<ArgumentsTypes>(args)...);
    } else {
      return ::mgis::internals::make_shared_impl<Type>(
          ctx, std::forward<ArgumentsTypes>(args)...);
    }
  }  // end of make_shared

  template <typename BaseType, typename Type, typename... ArgumentsTypes>
  std::shared_ptr<BaseType> make_shared_as(Context &ctx,
                                           ArgumentsTypes &&...args) noexcept
      requires std::is_base_of_v<BaseType, Type> &&
      ::mgis::internals::is_constructible<Type, ArgumentsTypes...> {
    if constexpr (::mgis::internals::is_constructible_with_context<
                      Type, ArgumentsTypes...>) {
      return ::mgis::internals::make_shared_as_impl<BaseType, Type>(
          ctx, ctx, std::forward<ArgumentsTypes>(args)...);
    } else {
      return ::mgis::internals::make_shared_as_impl<BaseType, Type>(
          ctx, std::forward<ArgumentsTypes>(args)...);
    }
  }  // end of make_shared_as

}  // end of namespace mgis

#endif /* LIB_MGIS_UTILITIES_CONSTRUCT_IXX */
