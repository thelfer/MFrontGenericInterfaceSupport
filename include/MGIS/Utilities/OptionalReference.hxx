/*!
 * \file   MGIS/Utilities/OptionalReference.hxx
 * \brief  This header declares the OptionalReference class
 * \author Thomas Helfer
 * \date   19/02/2026
 */

#ifndef LIB_MGIS_UTILITIES_OPTIONALREFERENCE_HXX
#define LIB_MGIS_UTILITIES_OPTIONALREFERENCE_HXX

#include <memory>
#include <cstddef>
#include <utility>
#include <cassert>
#include <type_traits>
#include "MGIS/Config.hxx"
#include "MGIS/InvalidResult.hxx"

namespace mgis {

  /*!
   * \brief a class that may contain a reference to a class
   *
   * \note this class was introduced by of the lack of optional
   * references prior to C++-26.
   * \note the design of the class is loosely modeled by the
   * std::experimental::observer_ptr proposal
   */
  template <typename ValueType>
  struct [[nodiscard]] OptionalReference {
    //
    using element_type = ValueType;
    using pointer = ValueType*;
    using reference = ValueType&;
    //
    constexpr OptionalReference() noexcept : ptr(nullptr) {}

    constexpr OptionalReference(std::nullptr_t) noexcept : ptr(nullptr) {}

    constexpr OptionalReference(pointer p) noexcept : ptr(p) {}

    template <typename ValueType2>
    requires(std::is_convertible<ValueType2*, ValueType*>::value)  //
        constexpr OptionalReference(
            OptionalReference<ValueType2> other) noexcept
        : ptr(other.get()) {}

    template <typename ValueType2>
    requires(std::is_convertible<ValueType2*, ValueType*>::value)  //
        constexpr OptionalReference(
            std::unique_ptr<ValueType2> const& other) noexcept
        : ptr(other.get()) {}

    template <typename ValueType2>
    requires(std::is_convertible<ValueType2*, ValueType*>::value)  //
        constexpr OptionalReference(
            std::shared_ptr<ValueType2> const& other) noexcept
        : ptr(other.get()) {}

    [[nodiscard]] constexpr pointer get() const noexcept { return this->ptr; }

    [[nodiscard]] constexpr reference operator*() const {
      assert(ptr != nullptr);
      return *ptr;
    }

    constexpr pointer operator->() const noexcept { return this->ptr; }

    [[nodiscard]] constexpr bool has_value() const noexcept {
      return this->ptr != nullptr;
    }

    [[nodiscard]] constexpr operator bool() const noexcept {
      return this->ptr != nullptr;
    }

    constexpr operator pointer() const noexcept { return this->ptr; }

    constexpr pointer release() noexcept {
      pointer p(ptr);
      this->reset();
      return p;
    }

    constexpr void reset(pointer p = nullptr) noexcept { this->ptr = p; }

    constexpr void swap(OptionalReference& other) noexcept {
      std::swap(ptr, other.ptr);
    }

   private:
    pointer ptr;
  };

  template <typename ValueType>
  void swap(OptionalReference<ValueType>& p1,
            OptionalReference<ValueType>& p2) noexcept {
    p1.swap(p2);
  }

  template <typename ValueType>
  [[nodiscard]] constexpr OptionalReference<ValueType> make_optional_reference(
      ValueType& p) noexcept {
    return OptionalReference<ValueType>(&p);
  }

  template <typename ValueType>
  [[nodiscard]] constexpr OptionalReference<ValueType> make_optional_reference(
      ValueType* p) noexcept {
    return OptionalReference<ValueType>(p);
  }

  template <typename ValueType1, typename ValueType2>
  [[nodiscard]] constexpr bool operator==(OptionalReference<ValueType1> p1,
                                          OptionalReference<ValueType2> p2) {
    return p1.get() == p2.get();
  }

  template <typename ValueType1, typename ValueType2>
  [[nodiscard]] constexpr bool operator!=(OptionalReference<ValueType1> p1,
                                          OptionalReference<ValueType2> p2) {
    return !(p1 == p2);
  }

  template <typename ValueType>
  [[nodiscard]] constexpr bool operator==(OptionalReference<ValueType> p,
                                          std::nullptr_t) noexcept {
    return !p;
  }

  template <typename ValueType>
  [[nodiscard]] constexpr bool operator==(
      std::nullptr_t, OptionalReference<ValueType> p) noexcept {
    return !p;
  }

  template <typename ValueType>
  [[nodiscard]] constexpr bool operator!=(OptionalReference<ValueType> p,
                                          std::nullptr_t) noexcept {
    return static_cast<bool>(p);
  }

  template <typename ValueType>
  [[nodiscard]] constexpr bool operator!=(
      std::nullptr_t, OptionalReference<ValueType> p) noexcept {
    return static_cast<bool>(p);
  }

  template <typename ValueType>
  [[nodiscard]] constexpr bool isInvalid(
      const OptionalReference<ValueType>& p) noexcept {
    return p == nullptr;
  }  // end of isInvalid

}  // end of namespace mgis

namespace std {

  template <typename ValueType>
  struct hash<::mgis::OptionalReference<ValueType>> {
    [[nodiscard]] size_t operator()(
        ::mgis::OptionalReference<ValueType> p) const noexcept {
      return hash<ValueType*>()(p.get());
    }
  };

}  // namespace std

namespace mgis::internal {

  //! \brief partial specialisation for optional references
  template <typename ValueType>
  struct InvalidValueTraits<mgis::OptionalReference<ValueType>> {
    static constexpr bool isSpecialized = true;
    static mgis::OptionalReference<ValueType> getValue() noexcept { return {}; }
  };

  //! \brief partial specialisation for optional reference
  template <typename T>
  struct OptionalTraits<OptionalReference<T>> {
    static constexpr bool isSpecialized = true;
    static constexpr T& deference(OptionalReference<T>&& v) noexcept {
      return *v;
    }
  };

}  // end of namespace mgis::internal

#endif /* LIB_MGIS_UTILITIES_OPTIONALREFERENCE_HXX */
