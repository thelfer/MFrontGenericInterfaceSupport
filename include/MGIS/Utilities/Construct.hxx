/*!
 * \file   MGIS/Utilities/Construct.hxx
 * \brief  This header declares utility functions to build objects handling any
 * exception thrown by a constructor \date   06/11/2022
 */

#ifndef LIB_MGIS_UTILITIES_CONSTRUCT_HXX
#define LIB_MGIS_UTILITIES_CONSTRUCT_HXX

#include <memory>
#include <optional>
#include <type_traits>
#include "MGIS/Context.hxx"
#include "MGIS/InvalidResult.hxx"

/*!
 * \def MGIS_CONSTRUCT(Type, e, ...)
 * \brief a simple wrapper around the `construct` function which handles the
 * declaration of the source location if needed
 *
 * \param[in] Type: type of the object built
 * \param[in, out] ctx: execution context
 */
#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
#define MGIS_CONSTRUCT(Type, e, ...) \
  ::mgis::construct<Type>(e, std::source_location::current(), __VA_ARGS__)
#else
#define MGIS_CONSTRUCT(Type, e, ...) ::mgis::construct<Type>(e, __VA_ARGS__)
#endif

/*!
 * \def MGIS_MAKE_UNIQUE(Type, e, ...)
 * \brief a simple wrapper around the `make_unique` function which handles the
 * declaration of the source location if needed.
 *
 * \param[in] Type: type of the object built
 * \param[in, out] ctx: execution context
 */
#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
#define MGIS_MAKE_UNIQUE(Type, e, ...) \
  ::mgis::make_unique<Type>(e, std::source_location::current(), __VA_ARGS__)
#else
#define MGIS_MAKE_UNIQUE(Type, e, ...) ::mgis::make_unique<Type>(e, __VA_ARGS__)
#endif

/*!
 * \def MGIS_MAKE_UNIQUE_AS(Type, e, ...)
 * \brief a simple wrapper around the `make_unique_as` function which handles
 * the declaration of the source location if needed
 *
 * \param[in] Type: type of the object built
 * \param[in, out] ctx: execution context
 */
#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
#define MGIS_MAKE_UNIQUE_AS(BaseType, Type, e, ...)                          \
  ::mgis::make_unique_as<BaseType, Type>(e, std::source_location::current(), \
                                         __VA_ARGS__)
#else
#define MGIS_MAKE_UNIQUE_AS(BaseType, Type, e, ...) \
  ::mgis::make_unique_as<BaseType, Type>(e, __VA_ARGS__)
#endif

/*!
 * \def MGIS_MAKE_SHARED(Type, e, ...)
 * \brief a simple wrapper around the `make_shared` function which handles the
 * declaration of the source location if needed
 *
 * \param[in] Type: type of the object built
 * \param[in, out] ctx: execution context
 */
#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
#define MGIS_MAKE_SHARED(Type, e, ...) \
  ::mgis::make_shared<Type>(e, std::source_location::current(), __VA_ARGS__)
#else
#define MGIS_MAKE_SHARED_AS(BaseType, Type, e, ...) \
  ::mgis::make_shared_as<BaseType, Type>(e, __VA_ARGS__)
#endif

/*!
 * \def MGIS_MAKE_SHARED_AS(Type, e, ...)
 * \brief a simple wrapper around the `make_shared_as` function which handles
 * the declaration of the source location if needed.
 *
 * \param[in] Type: type of the object built
 * \param[in, out] ctx: execution context
 */
#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
#define MGIS_MAKE_SHARED_AS(BaseType, Type, e, ...)                          \
  ::mgis::make_shared_as<BaseType, Type>(e, std::source_location::current(), \
                                         __VA_ARGS__)
#else
#define MGIS_MAKE_SHARED(Type, e, ...) ::mgis::make_shared<Type>(e, __VA_ARGS__)
#endif

/*!
 * \brief an helper macro to build a variable the constructor of which may throw
 *
 * \param[in] Type: type of the object built
 * \param[in] v: name of the variable built
 * \param[in, out] ctx: execution context
 */
#define MGIS_TRY_CONSTRUCT(Type, v, e, ...)        \
  std::optional<Type> MGIS_TEMPORARY_VARIABLE(v) = \
      MGIS_CONSTRUCT(Type, e, __VA_ARGS__);        \
  if (!MGIS_TEMPORARY_VARIABLE(v).has_value()) {   \
    return InvalidResult{};                        \
  }                                                \
  Type &v = *(MGIS_TEMPORARY_VARIABLE(v));

/*!
 * \brief an helper macro to build an unique pointer for type the constructor of
 * which may throw
 *
 * \param[in] Type: type of the object built
 * \param[in] v: name of the variable built
 * \param[in, out] ctx: execution context
 */
#define MGIS_TRY_MAKE_UNIQUE(Type, v, e, ...)                       \
  std::unique_ptr<Type> v = MGIS_MAKE_UNIQUE(Type, e, __VA_ARGS__); \
  if (v.get() == nullptr) {                                         \
    return InvalidResult{};                                         \
  }

/*!
 * \brief an helper macro to build an unique pointer for type the constructor of
 * which may throw
 *
 * \param[in] BaseType: pointer to a base type
 * \param[in] Type: type of the object built
 * \param[in] v: name of the variable built
 * \param[in, out] ctx: execution context
 */
#define MGIS_TRY_MAKE_UNIQUE_AS(BaseType, Type, v, e, ...) \
  std::unique_ptr<BaseType> v =                            \
      MGIS_MAKE_UNIQUE_AS(BaseType, Type, e, __VA_ARGS__); \
  if (v.get() == nullptr) {                                \
    return InvalidResult{};                                \
  }

/*!
 * \brief an helper macro to build an shared pointer for type the constructor of
 * which may throw
 *
 * \param[in] Type: type of the object built
 * \param[in] v: name of the variable built
 * \param[in, out] ctx: execution context
 */
#define MGIS_TRY_MAKE_SHARED(Type, v, e, ...)                       \
  std::shared_ptr<Type> v = MGIS_MAKE_SHARED(Type, e, __VA_ARGS__); \
  if (v.get() == nullptr) {                                         \
    return InvalidResult{};                                         \
  }

/*!
 * \brief an helper macro to build an shared pointer for type the constructor of
 * which may throw
 *
 * \param[in] BaseType: pointer to a base type
 * \param[in] Type: type of the object built
 * \param[in] v: name of the variable built
 * \param[in, out] ctx: execution context
 */
#define MGIS_TRY_MAKE_SHARED_AS(BaseType, Type, v, e, ...) \
  std::shared_ptr<BaseType> v =                            \
      MGIS_MAKE_SHARED_AS(BaseType, Type, e, __VA_ARGS__); \
  if (v.get() == nullptr) {                                \
    return InvalidResult{};                                \
  }

namespace mgis {

#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION

  /*!
   * \brief try to build an object of the given type.
   *
   * The constructor of the object may throw. In this case, the exception thrown
   * is catched and the information hold by the exception is added to the list
   * of registered error messages.
   *
   * \tparam Type: type to build.
   * \tparam ArgumentsTypes: types of the arguments passed to the constructor.
   * \param[in, out] ctx: execution context
   * \param[in] args: arguments passed to the constructor
   */
  template <typename Type, typename... ArgumentsTypes>
  [[nodiscard]] std::optional<Type> construct(
      Context &,
      const std::source_location &,
      ArgumentsTypes &&...) noexcept requires
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...>;
  /*!
   * \brief try to build a unique pointer holding the given type.
   *
   * The constructor of the object may throw. In this case, the exception thrown
   * is catched and the information hold by the exception is added to the list
   * of registered error messages.
   *
   * \tparam BaseType: base type used to store the build object.
   * \tparam Type: type to build.
   * \tparam ArgumentsTypes: types of the arguments passed to the constructor.
   * \param[in, out] ctx: execution context
   * \param[in] args: arguments passed to the constructor
   */
  template <typename Type, typename... ArgumentsTypes>
  [[nodiscard]] std::unique_ptr<Type> make_unique(
      Context &,
      const std::source_location &,
      ArgumentsTypes &&...) noexcept requires
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...>;
  /*!
   * \brief try to build a unique pointer of a base type holding the given type.
   *
   * The constructor of the object may throw. In this case, the exception thrown
   * is catched and the information hold by the exception is added to the list
   * of registered error messages.
   *
   * \tparam BaseType: base type used to store the build object.
   * \tparam Type: type to build.
   * \tparam ArgumentsTypes: types of the arguments passed to the constructor.
   * \param[in, out] ctx: execution context
   * \param[in] args: arguments passed to the constructor
   */
  template <typename BaseType, typename Type, typename... ArgumentsTypes>
  [[nodiscard]] std::unique_ptr<BaseType> make_unique_as(
      Context &, const std::source_location &, ArgumentsTypes &&...) noexcept
      requires std::is_base_of_v<BaseType, Type> &&
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...>;
  /*!
   * \brief try to build a shared pointer holding the given type.
   *
   * The constructor of the object may throw. In this case, the exception thrown
   * is catched and the information hold by the exception is added to the list
   * of registered error messages.
   *
   * \tparam Type: type to build.
   * \tparam ArgumentsTypes: types of the arguments passed to the constructor.
   * \param[in, out] ctx: execution context
   * \param[in] args: arguments passed to the constructor
   */
  template <typename Type, typename... ArgumentsTypes>
  [[nodiscard]] std::shared_ptr<Type> make_shared(
      Context &,
      const std::source_location &,
      ArgumentsTypes &&...) noexcept requires
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...>;
  /*!
   * \brief try to build a shared pointer of a base type holding the given type.
   *
   * The constructor of the object may throw. In this case, the exception thrown
   * is catched and the information hold by the exception is added to the list
   * of registered error messages.
   *
   * \tparam BaseType: base type used to store the build object.
   * \tparam Type: type to build.
   * \tparam ArgumentsTypes: types of the arguments passed to the constructor.
   * \param[in, out] ctx: execution context
   * \param[in] args: arguments passed to the constructor
   */
  template <typename BaseType, typename Type, typename... ArgumentsTypes>
  [[nodiscard]] std::shared_ptr<BaseType> make_shared_as(
      Context &, const std::source_location &, ArgumentsTypes &&...) noexcept
      requires std::is_base_of_v<BaseType, Type> &&
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...>;

#endif

  /*!
   * \brief try to build an object of the given type.
   *
   * The constructor of the object may throw. In this case, the exception thrown
   * is catched and the information hold by the exception is added to the list
   * of registered error messages.
   *
   * \tparam Type: type to build.
   * \tparam ArgumentsTypes: types of the arguments passed to the constructor.
   * \param[in, out] ctx: execution context
   * \param[in] args: arguments passed to the constructor
   */
  template <typename Type, typename... ArgumentsTypes>
  [[nodiscard]] std::optional<Type> construct(
      Context &, ArgumentsTypes &&...) noexcept requires
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...>;
  /*!
   * \brief try to build a unique pointer holding the given type.
   *
   * The constructor of the object may throw. In this case, the exception thrown
   * is catched and the information hold by the exception is added to the list
   * of registered error messages.
   *
   * \tparam BaseType: base type used to store the build object.
   * \tparam Type: type to build.
   * \tparam ArgumentsTypes: types of the arguments passed to the constructor.
   * \param[in, out] ctx: execution context
   * \param[in] args: arguments passed to the constructor
   */
  template <typename Type, typename... ArgumentsTypes>
  [[nodiscard]] std::unique_ptr<Type> make_unique(
      Context &, ArgumentsTypes &&...) noexcept requires
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...>;
  /*!
   * \brief try to build a unique pointer of a base type holding the given type.
   *
   * The constructor of the object may throw. In this case, the exception thrown
   * is catched and the information hold by the exception is added to the list
   * of registered error messages.
   *
   * \tparam BaseType: base type used to store the build object.
   * \tparam Type: type to build.
   * \tparam ArgumentsTypes: types of the arguments passed to the constructor.
   * \param[in, out] ctx: execution context
   * \param[in] args: arguments passed to the constructor
   */
  template <typename BaseType, typename Type, typename... ArgumentsTypes>
  [[nodiscard]] std::unique_ptr<BaseType> make_unique_as(
      Context &, ArgumentsTypes &&...) noexcept requires
      std::is_base_of_v<BaseType, Type> &&
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...>;
  /*!
   * \brief try to build a shared pointer holding the given type.
   *
   * The constructor of the object may throw. In this case, the exception thrown
   * is catched and the information hold by the exception is added to the list
   * of registered error messages.
   *
   * \tparam Type: type to build.
   * \tparam ArgumentsTypes: types of the arguments passed to the constructor.
   * \param[in, out] ctx: execution context
   * \param[in] args: arguments passed to the constructor
   */
  template <typename Type, typename... ArgumentsTypes>
  [[nodiscard]] std::shared_ptr<Type> make_shared(
      Context &, ArgumentsTypes &&...) noexcept requires
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...>;
  /*!
   * \brief try to build a shared pointer of a base type holding the given type.
   *
   * The constructor of the object may throw. In this case, the exception thrown
   * is catched and the information hold by the exception is added to the list
   * of registered error messages.
   *
   * \tparam BaseType: base type used to store the build object.
   * \tparam Type: type to build.
   * \tparam ArgumentsTypes: types of the arguments passed to the constructor.
   * \param[in, out] ctx: execution context
   * \param[in] args: arguments passed to the constructor
   */
  template <typename BaseType, typename Type, typename... ArgumentsTypes>
  [[nodiscard]] std::shared_ptr<BaseType> make_shared_as(
      Context &, ArgumentsTypes &&...) noexcept requires
      std::is_base_of_v<BaseType, Type> &&
      std::is_constructible_v<std::remove_const_t<Type>, ArgumentsTypes...>;

}  // end of namespace mgis

#include "MGIS/Utilities/Construct.ixx"

#endif /* LIB_MGIS_UTILITIES_CONSTRUCT_HXX */
