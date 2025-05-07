/*!
 * \file   MGIS/InvalidResult.hxx
 * \brief  This file declares the `InvalidResult` class.
 * \date   15/11/2022
 */

#ifndef LIB_MGIS_INVALIDRESULT_HXX
#define LIB_MGIS_INVALIDRESULT_HXX 1

#include <memory>
#include <optional>

namespace mgis::internal {

/*!
 * \brief metafunction returning an invalid value associated with a type
 *
 * The InvalidResult class can be converted to any type for which this
 * class is specialized.
 *
 * Specialisation of this class must expose:
 *
 * - a public static constexpr data member named `isSpecialized` egal to true.
 * - a public static constexpr method named `getValue` returning the invalid value.
 */
template<typename Type>
struct InvalidValueTraits {
  //
  static constexpr bool isSpecialized = false;
};

}    // end of namespace mgis::internal

namespace mgis {

/*!
 * \brief a class convertible to many values used as invalid result.
 *
 * This class has been designed to make the `MANTA_TRY*` macros
 * compatible with any function/methods respecting the guideline
 * of error handling.
 *
 * Internally, these macros can return an object of the
 * `InvalidResult` type when an error is detected, which will be
 * automatically converted to the return type of the current
 * function/method.
 */
struct InvalidResult {
  /*!
   * \brief conversion to the invalid value to any type for which
   * a specialisation of the InvalidValue class exists
   */
  template <typename Type>
  [[nodiscard]] constexpr operator Type() &&noexcept
      requires(internal::InvalidValueTraits<Type>::isSpecialized);
  //! \brief implicit conversion to std::optional<bool>
  [[nodiscard]] constexpr operator std::optional<bool>() && noexcept;

};    // end of class InvalidResult
/*!
 * \return if the given boolean is true
 * \param[in] b: value
 */
[[nodiscard]] constexpr bool isInvalid(const bool) noexcept;
/*!
 * \return if the given optional is invalid (does not contain an object)
 * \param[in] o: optional value
 */
template <typename T>
[[nodiscard]] constexpr bool isInvalid(const std::optional<T> &) noexcept;
/*!
 * \return if the given pointer is null (i.e. is empty)
 * \param[in] p: pointer
 */
template <typename T>
[[nodiscard]] constexpr bool isInvalid(T *const) noexcept;
/*!
 * \return if the given std::unique_ptr is invalid (i.e. is empty)
 * \param[in] p: pointer
 */
template<typename T>
[[nodiscard]] constexpr bool isInvalid(const std::unique_ptr<T> &) noexcept;
/*!
 * \return if the given std::shared_ptr is invalid (i.e. is empty)
 * \param[in] p: pointer
 */
template<typename T>
[[nodiscard]] constexpr bool isInvalid(const std::shared_ptr<T> &) noexcept;
/*!
 * \return if one of the given objects is invalid
 * \param[in] o: objects tested
 */
[[nodiscard]] constexpr bool areInvalid(const auto &...) noexcept;
/*!
 * \return if the given object is valid
 * \param[in] o: object tested
 */
[[nodiscard]] constexpr bool isValid(const auto &) noexcept;
/*!
 * \return if all the given objects are valid
 * \param[in] o: objects tested
 */
[[nodiscard]] constexpr bool areValid(const auto &...) noexcept;

}    // end of namespace mgis

#include "MGIS/InvalidResult.ixx"

#endif /* LIB_MGIS_INVALIDRESULT_HXX */
