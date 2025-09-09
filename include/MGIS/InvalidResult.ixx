/*!
 * \file   MGIS/InvalidResult.ixx
 * \brief  This file implements the inline methods of the `InvalidResult` class.
 * \date   15/11/2022
 */

#ifndef LIB_MGIS_INVALIDRESULT_IXX
#define LIB_MGIS_INVALIDRESULT_IXX 1

namespace mgis::internal {

  //! \brief partial specialisation for boolean values
  template <>
  struct InvalidValueTraits<bool> {
    static constexpr bool isSpecialized = true;
    static bool getValue() noexcept { return false; }
  };

  //! \brief partial specialisation for std::optional
  template <typename T>
  struct InvalidValueTraits<std::optional<T>> {
    static constexpr bool isSpecialized = true;
    static constexpr auto getValue() noexcept { return std::optional<T>{}; }
  };

  //! \brief partial specialisation for std::optional
  template <typename T>
  struct InvalidValueTraits<std::optional<const T>> {
    static constexpr bool isSpecialized = true;
    static constexpr auto getValue() noexcept {
      return std::optional<const T>{};
    }
  };

  //! \brief partial specialisation for std::optional<bool>
  template <>
  struct InvalidValueTraits<std::optional<bool>> {
    static constexpr bool isSpecialized = true;
    static constexpr auto getValue() noexcept { return std::optional<bool>{}; }
  };

  //! \brief partial specialisation for std::unique_ptr
  template <typename T>
  struct InvalidValueTraits<std::unique_ptr<T>> {
    static constexpr bool isSpecialized = true;
    static constexpr auto getValue() noexcept { return std::unique_ptr<T>{}; }
  };

  //! \brief partial specialisation for std::shared_ptr
  template <typename T>
  struct InvalidValueTraits<std::shared_ptr<T>> {
    static constexpr bool isSpecialized = true;
    static constexpr auto getValue() noexcept { return std::shared_ptr<T>{}; }
  };

}  // end of namespace mgis::internal

namespace mgis {

  template <typename Type>
  [[nodiscard]] constexpr InvalidResult::operator Type() &&noexcept
      requires(internal::InvalidValueTraits<Type>::isSpecialized) {
    return internal::InvalidValueTraits<Type>::getValue();
  }  // end of operator Type

  [[nodiscard]] constexpr InvalidResult::operator std::optional<bool>()
      &&noexcept {
    return {};
  }  // end of operator std::optional<bool>

  constexpr bool isInvalid(const bool b) noexcept {
    return !b;
  }  // end of is Invalid

  template <typename T>
  constexpr bool isInvalid(T *const ptr) noexcept {
    return ptr == nullptr;
  }  // end of is Invalid

  template <typename T>
  constexpr bool isInvalid(const std::optional<T> &o) noexcept {
    return !o.has_value();
  }  // end of is Invalid

  template <typename T>
  constexpr bool isInvalid(const std::unique_ptr<T> &p) noexcept {
    return p.get() == nullptr;
  }  // end of is Invalid

  template <typename T>
  constexpr bool isInvalid(const std::shared_ptr<T> &p) noexcept {
    return p.get() == nullptr;
  }  // end of is isInvalid

  constexpr bool areInvalid(const auto &...o) noexcept {
    return (... || isInvalid(o));
  }  // end of is isInvalid

  constexpr bool isValid(const auto &o) noexcept {
    return !isInvalid(o);
  }  // end of is isValid

  constexpr bool areValid(const auto &...o) noexcept {
    return (... && isValid(o));
  }  // end of isValid

}  // end of namespace mgis

#endif /* LIB_MGIS_INVALIDRESULT_IXX */
