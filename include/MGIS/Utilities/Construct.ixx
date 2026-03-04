/*!
 * \file   MIGS/Utilities/Construct.ixx
 * \brief  This file implements the template functions declared
 *         in the `MGIS/Construct.hx` header.
 * \date   04/11/2022
 */

#ifndef LIB_MGIS_UTILITIES_CONSTRUCT_IXX
#define LIB_MGIS_UTILITIES_CONSTRUCT_IXX

namespace mgis {

#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION

  template <typename Type, typename... ArgumentsTypes>
  std::optional<Type> construct(Context &ctx,
                                const std::source_location &l,
                                ArgumentsTypes &&...args) noexcept requires
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...> {
    using NonConstType = std::remove_const_t<Type>;
    if constexpr (std::is_nothrow_constructible_v<NonConstType,
                                                  ArgumentsTypes...>) {
      return makestd::optional<Type>(std::forward<ArgumentsTypes>(args)...);
    } else {
      try {
        return makestd::optional<Type>(std::forward<ArgumentsTypes>(args)...);
      } catch (...) {
        registerExceptionInErrorBacktrace(ctx, l);
      }
      return {};
    }
  }  // end of construct

  template <typename Type, typename... ArgumentsTypes>
  std::unique_ptr<Type> make_unique(Context &ctx,
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
  std::unique_ptr<BaseType> make_unique_as(Context &ctx,
                                           const std::source_location &l,
                                           ArgumentsTypes &&...args) noexcept
      requires std::is_base_of_v<BaseType, Type> &&
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
  std::shared_ptr<Type> make_shared(Context &ctx,
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
  std::shared_ptr<BaseType> make_shared_as(Context &ctx,
                                           const std::source_location &l,
                                           ArgumentsTypes &&...args) noexcept
      requires std::is_base_of_v<BaseType, Type> &&
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
  std::optional<Type> construct(Context &ctx,
                                ArgumentsTypes &&...args) noexcept requires
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...> {
    using NonConstType = std::remove_const_t<Type>;
    if constexpr (std::is_nothrow_constructible_v<NonConstType,
                                                  ArgumentsTypes...>) {
      return Type{std::forward<ArgumentsTypes>(args)...};
    } else {
      try {
        return std::optional<Type>(std::forward<ArgumentsTypes>(args)...);
      } catch (...) {
        std::ignore =
            registerExceptionInErrorBacktraceWithoutSourceLocation(ctx);
      }
      return {};
    }
  }  // end of construct

  template <typename Type, typename... ArgumentsTypes>
  std::unique_ptr<Type> make_unique(Context &ctx,
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
        std::ignore =
            registerExceptionInErrorBacktraceWithoutSourceLocation(ctx);
      }
      return {};
    }
  }  // end of make_unique

  template <typename BaseType, typename Type, typename... ArgumentsTypes>
  std::unique_ptr<BaseType> make_unique_as(Context &ctx,
                                           ArgumentsTypes &&...args) noexcept
      requires std::is_base_of_v<BaseType, Type> &&
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
  std::shared_ptr<Type> make_shared(Context &ctx,
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
        std::ignore =
            registerExceptionInErrorBacktraceWithoutSourceLocation(ctx);
      }
      return {};
    }
  }  // end of make_shared

  template <typename BaseType, typename Type, typename... ArgumentsTypes>
  std::shared_ptr<BaseType> make_shared_as(Context &ctx,
                                           ArgumentsTypes &&...args) noexcept
      requires std::is_base_of_v<BaseType, Type> &&
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

}  // namespace mgis

#endif /* LIB_MGIS_UTILITIES_CONSTRUCT_IXX */