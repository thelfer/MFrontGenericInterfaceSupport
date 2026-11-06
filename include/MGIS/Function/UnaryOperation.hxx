/*!
 * \file   MGIS/Function/UnaryOperation.hxx
 * \brief  This file declares the `UnaryOperation` and `UnaryOperation2`
 * classes, which build a new evaluator by applying a unary operation to the
 * values returned by an existing evaluator, the `unary_operation` and
 * `transform` functions used to build such evaluators, the customization
 * points used to implement unary operations on the supported result types and
 * the predefined modifiers `absolute_value`, `maximum_component`,
 * `minimum_component`, `negate`, `multiply_by_scalar` and
 * `divide_by_scalar`.
 * \author Thomas Helfer
 * \date   09/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_UNARYOPERATION_HXX
#define LIB_MGIS_FUNCTION_UNARYOPERATION_HXX

#include <span>
#include <algorithm>
#include <type_traits>

#ifdef MGIS_HAVE_TFEL
#include "TFEL/Math/General/BasicOperations.hxx"
#include "TFEL/Math/General/UnaryResultType.hxx"
#include "TFEL/Math/General/ResultType.hxx"
#endif /* MGIS_HAVE_TFEL */

#include "MGIS/Function/EvaluatorModifierBase.hxx"

namespace mgis::function {

  namespace internals {

    /*!
     * \brief helper base class giving the number of components of the values
     * returned by a unary operation applied to an evaluator
     *
     * The number of components is the compile-time size of the type returned
     * by the callable when it is called with the values returned by the
     * evaluator. A static assertion is raised when this number is only known
     * at run time, as evaluators are required to know it at compile time.
     *
     * \tparam CallableType: type of the unary operation
     * \tparam EvaluatorType: type of the evaluator whose values are modified
     */
    template <typename CallableType, EvaluatorConcept EvaluatorType>
    struct UnaryOperationBase {
      /*!
       * \brief return the number of components
       * \return number of components of the values returned by the unary
       * operation, which is known at compile time
       */
      static constexpr size_type getNumberOfComponents() MGIS_NOEXCEPT {
        using UnaryOperationResult =
            std::invoke_result_t<CallableType, evaluator_result<EvaluatorType>>;
        return compile_time_size<UnaryOperationResult>;
      }
      static_assert(getNumberOfComponents() != dynamic_extent);
    };

  }  // end of namespace internals

  /*!
   * \brief an evaluator applying a unary operation to the values returned by
   * another evaluator
   *
   * Such an evaluator is usually built with the `unary_operation` or
   * `transform` functions, which can either be applied directly to an
   * evaluator or be composed with other modifiers using `operator|`:
   *
   * \code{.cpp}
   * const auto e2 = e | transform([](const real v) { return v * v; });
   * \endcode
   *
   * The callable object is stored in the instance. The number of components
   * of the returned values is the compile-time size of the type returned by
   * the callable, which must thus be known at compile time. Only the call
   * operators which match the way the underlying evaluator can be called are
   * available, as they are inherited from `EvaluatorModifierBase`.
   *
   * \tparam CallableType: type of the unary operation, which must be copy
   * constructible and callable with the values returned by the evaluator
   * \tparam EvaluatorType: type of the evaluator whose values are modified
   */
  template <typename CallableType, EvaluatorConcept EvaluatorType>
  requires((std::is_copy_constructible_v<CallableType>)&&(
      std::invocable<CallableType,
                     evaluator_result<EvaluatorType>>))  //
      struct UnaryOperation
      : EvaluatorModifierBase<UnaryOperation<CallableType, EvaluatorType>,
                              EvaluatorType>,
        internals::UnaryOperationBase<CallableType, EvaluatorType> {
    /*!
     * \brief constructor
     * \param[in] c: unary operation, applied to the values returned by the
     * evaluator, a copy of it is stored in the instance
     * \param[in] e: evaluator providing those values, a copy of it is stored
     * in the instance
     */
    MGIS_HOST_DEVICE constexpr UnaryOperation(const CallableType&,
                                              const EvaluatorType&);
    //
    using internals::UnaryOperationBase<CallableType,
                                        EvaluatorType>::getNumberOfComponents;
    /*!
     * \brief apply the modifier
     * \param[in] values: values computed by the underlying evaluator at a
     * given integration point
     * \return values computed by calling the stored callable with the values
     * computed by the underlying evaluator at a given integration point
     * \note this function is called by the call operators inherited from
     * `EvaluatorModifierBase`
     */
    MGIS_HOST_DEVICE constexpr auto apply(
        const evaluator_result<EvaluatorType>&) const;

   private:
    //! \brief unary operation applied to the values of the evaluator
    CallableType modifier;
  };

  /*!
   * \brief an evaluator applying a unary operation to the values returned by
   * another evaluator, without storing this operation
   *
   * This class behaves as `UnaryOperation` except that the callable object is
   * not stored in the instance: it must be trivially default constructible and
   * a default constructed object is created each time the values are computed.
   * Such a modifier is thus as light as the underlying evaluator, which allows
   * `constexpr` modifiers such as `absolute_value` or `negate` to be defined.
   *
   * The constructor inherited from `EvaluatorModifierBase` only takes the
   * underlying evaluator as an argument.
   *
   * \tparam CallableType: type of the unary operation, which must be trivially
   * default constructible and callable with the values returned by the
   * evaluator
   * \tparam EvaluatorType: type of the evaluator whose values are modified
   */
  template <typename CallableType, EvaluatorConcept EvaluatorType>
  requires((std::is_trivially_default_constructible_v<CallableType>)&&(
      std::invocable<CallableType,
                     evaluator_result<EvaluatorType>>))  //
      struct UnaryOperation2
      : EvaluatorModifierBase<UnaryOperation2<CallableType, EvaluatorType>,
                              EvaluatorType>,
        internals::UnaryOperationBase<CallableType, EvaluatorType> {
    // inheriting constructor
    using EvaluatorModifierBase<UnaryOperation2,
                                EvaluatorType>::EvaluatorModifierBase;
    //
    using internals::UnaryOperationBase<CallableType,
                                        EvaluatorType>::getNumberOfComponents;
    /*!
     * \brief apply the modifier
     * \param[in] values: values computed by the underlying evaluator at a
     * given integration point
     * \return values computed by calling a default constructed object of
     * type `CallableType` with the values computed by the underlying
     * evaluator at a given integration point
     * \note this function is called by the call operators inherited from
     * `EvaluatorModifierBase`
     */
    MGIS_HOST_DEVICE constexpr auto apply(
        const evaluator_result<EvaluatorType>&) const;
  };

  /*!
   * \brief return the number of components of the values returned by a unary
   * operation
   * \tparam CallableType: type of the unary operation
   * \tparam EvaluatorType: type of the evaluator whose values are modified
   * \param[in] e: unary operation
   * \return number of components of the values returned by the operation
   */
  template <typename CallableType, EvaluatorConcept EvaluatorType>
  constexpr mgis::size_type getNumberOfComponents(
      const UnaryOperation<CallableType, EvaluatorType>&) MGIS_NOEXCEPT;
  /*!
   * \brief return the number of components of the values returned by a unary
   * operation which does not store its callable object
   * \tparam CallableType: type of the unary operation
   * \tparam EvaluatorType: type of the evaluator whose values are modified
   * \param[in] e: unary operation
   * \return number of components of the values returned by the operation
   */
  template <typename CallableType, EvaluatorConcept EvaluatorType>
  constexpr mgis::size_type getNumberOfComponents(
      const UnaryOperation2<CallableType, EvaluatorType>&) MGIS_NOEXCEPT;

  namespace internals {

    /*!
     * \brief a modifier building a `UnaryOperation` from a callable object
     * which is stored in the modifier itself
     *
     * \tparam CallableType: type of the unary operation, which must be copy
     * constructible
     */
    template <typename CallableType>
    struct unary_operation_modifier {
      //! \brief this alias allows to match the Evaluator Modifier concept
      using Tag = ::mgis::function::EvaluatorModifierTag;
      /*!
       * \brief constructor
       * \param[in] c: unary operation, a copy of it is stored in the modifier
       */
      constexpr unary_operation_modifier(const CallableType&);
      /*!
       * \brief apply the unary operation to an evaluator
       * \tparam EvaluatorType: type of the evaluator, which may be a qualified
       * type
       * \param[in] e: evaluator whose values are modified
       * \return a `UnaryOperation` storing a copy of the callable object and
       * the given evaluator
       */
      template <typename EvaluatorType>
      constexpr auto operator()(EvaluatorType&&) const
          requires((EvaluatorConcept<std::decay_t<EvaluatorType>>)&&(
              std::invocable<CallableType,
                             evaluator_result<std::decay_t<EvaluatorType>>>));

     private:
      //! \brief unary operation applied to the values of the evaluator
      CallableType modifier;
    };

    /*!
     * \brief a modifier building a `UnaryOperation2` from a callable object
     * which is not stored in the modifier
     *
     * This helper is a stateless type, which allows the predefined modifiers
     * of this file to be defined as `constexpr` global variables.
     *
     * \tparam CallableType: type of the unary operation, which must be
     * trivially default constructible
     */
    template <typename CallableType>
    struct unary_operation_modifier2_impl {
      //! \brief this alias allows to match the Evaluator Modifier concept
      using Tag = ::mgis::function::EvaluatorModifierTag;
      //
      /*!
       * \brief apply the unary operation to an evaluator
       * \tparam EvaluatorType: type of the evaluator, which may be a qualified
       * type
       * \param[in] e: evaluator whose values are modified
       * \return a `UnaryOperation2` wrapping the given evaluator, the callable
       * object being created when the values are computed
       */
      template <typename EvaluatorType>
      MGIS_HOST_DEVICE constexpr auto operator()(EvaluatorType&&) const
          requires((EvaluatorConcept<std::decay_t<EvaluatorType>>)&&(
              std::invocable<CallableType,
                             evaluator_result<std::decay_t<EvaluatorType>>>));
    };

    /*!
     * \brief create a modifier from a callable object which is only used to
     * deduce the type of the unary operation
     *
     * \tparam CallableType: type of the unary operation, which must be
     * trivially default constructible
     * \param[in] c: callable object which is discarded, as the operation is
     * recreated when the values are computed
     * \return a stateless modifier of type
     * `unary_operation_modifier2_impl<CallableType>`
     */
    template <typename CallableType>
    MGIS_HOST_DEVICE constexpr auto unary_operation_modifier2(CallableType);

  }  // namespace internals

  /*!
   * \brief create a modifier applying a unary operation to the values returned
   * by an evaluator
   *
   * The returned modifier can be applied to an evaluator using `operator|` or
   * composed with other modifiers, as in
   * `e | unary_operation(c1) | unary_operation(c2)`.
   *
   * \tparam CallableType: type of the unary operation, which may be a qualified
   * type
   * \param[in] c: unary operation, a copy of it is stored in the modifier
   * \return a `internals::unary_operation_modifier` object
   */
  template <typename CallableType>
  constexpr auto unary_operation(CallableType&&);

  /*!
   * \brief apply a unary operation to the values returned by an evaluator
   *
   * \tparam CallableType: type of the unary operation, which may be a qualified
   * type
   * \tparam EvaluatorType: type of the evaluator
   * \param[in] c: unary operation, which must be callable with the values
   * returned by the evaluator
   * \param[in] e: evaluator whose values are modified
   * \return a `UnaryOperation` object, i.e. an evaluator returning the values
   * computed by the callable object
   */
  template <typename CallableType, EvaluatorConcept EvaluatorType>
  constexpr auto unary_operation(CallableType&&, const EvaluatorType&);

  /*!
   * \brief create a modifier transforming the values returned by an evaluator
   *
   * This function is an alias of the `unary_operation` function, which is often
   * more readable when the operation is a user defined lambda:
   *
   * \code{.cpp}
   * const auto e2 = e | transform([](const real v) { return v * v; });
   * \endcode
   *
   * \tparam CallableType: type of the unary operation, which may be a qualified
   * type
   * \param[in] c: unary operation, a copy of it is stored in the modifier
   * \return a `internals::unary_operation_modifier` object
   */
  template <typename CallableType>
  constexpr auto transform(CallableType&&);

  /*!
   * \brief transform the values returned by an evaluator
   *
   * This function is an alias of the `unary_operation` function taking an
   * evaluator as second argument.
   *
   * \tparam CallableType: type of the unary operation, which may be a qualified
   * type
   * \tparam EvaluatorType: type of the evaluator
   * \param[in] c: unary operation, which must be callable with the values
   * returned by the evaluator
   * \param[in] e: evaluator whose values are modified
   * \return a `UnaryOperation` object, i.e. an evaluator returning the values
   * computed by the callable object
   */
  template <typename CallableType, EvaluatorConcept EvaluatorType>
  constexpr auto transform(CallableType&&, const EvaluatorType&);

}  // end of namespace mgis::function

#include "MGIS/Function/UnaryOperation.ixx"

namespace mgis::function::customization_points {

  /*!
   * \brief customization point giving the absolute value of a value which can
   * be returned by an evaluator
   *
   * The primary template is only declared. Specializations are provided for the
   * result types supported by the library (`real`, fixed size constant spans
   * and arrays of scalar values) and new specializations can be added to
   * support other result types.
   *
   * \tparam T: type of the value whose absolute value is requested
   */
  template <typename T>
  struct AbsoluteValue;

  //! \brief specialization for scalar values
  template <>
  struct AbsoluteValue<real> {
    /*!
     * \brief compute the absolute value of a scalar value
     * \param[in] v: scalar value
     * \return `v` if `v` is positive or null, `-v` otherwise
     */
    [[nodiscard]] static constexpr real exe(const real& v) MGIS_NOEXCEPT {
      return v < 0 ? -v : v;
    }
  };

  /*!
   * \brief specialization for constant fixed sized spans
   *
   * The absolute values cannot be returned as a view of the given span since
   * the values it points to are constant: they are copied instead.
   *
   * \tparam N: number of components of the span, which must be strictly
   * positive and known at compile time
   */
  template <std::size_t N>
  requires((N > 0) && (N != std::dynamic_extent))  //
      struct AbsoluteValue<std::span<const real, N>> {
    /*!
     * \brief compute the absolute value of each component of a span
     * \param[in] v: span holding `N` scalar values
     * \return a `std::array<real, N>` containing the absolute value of each
     * component of the given span
     */
    [[nodiscard]] static constexpr auto exe(const std::span<const real, N>& v)
        MGIS_NOEXCEPT {
      auto r = std::array<real, N>{};
      for (std::size_t i = 0; i != N; ++i) {
        r[i] = std::abs(v[i]);
      }
      return r;
    }
  };

  /*!
   * \brief specialization for fixed sized arrays of scalar values
   *
   * The values are copied in the returned array.
   *
   * \tparam N: number of components of the array, which must be strictly
   * positive and known at compile time
   */
  template <std::size_t N>
  requires((N > 0) && (N != std::dynamic_extent))  //
      struct AbsoluteValue<std::array<real, N>> {
    /*!
     * \brief compute the absolute value of each component of an array
     * \param[in] v: array holding `N` scalar values
     * \return a `std::array<real, N>` containing the absolute value of each
     * component of the given array
     */
    [[nodiscard]] static constexpr auto exe(const std::array<real, N>& v)
        MGIS_NOEXCEPT {
      auto r = std::array<real, N>{};
      for (std::size_t i = 0; i != N; ++i) {
        r[i] = std::abs(v[i]);
      }
      return r;
    }
  };

  /*!
   * \brief customization point giving the largest component of a value which
   * can be returned by an evaluator
   *
   * The primary template is only declared. Specializations are provided for the
   * result types supported by the library (`real`, fixed size constant spans
   * and arrays of scalar values) and new specializations can be added to
   * support other result types.
   *
   * \tparam T: type of the value whose largest component is requested
   */
  template <typename T>
  struct MaximumComponent;

  /*!
   * \brief specialization for scalar values
   *
   * A scalar value has a single component: the largest component is the value
   * itself.
   */
  template <>
  struct MaximumComponent<real> {
    /*!
     * \brief return the largest component of a scalar value
     * \param[in] v: scalar value
     * \return the given value
     */
    [[nodiscard]] static constexpr real exe(const real& v) MGIS_NOEXCEPT {
      return v;
    }
  };

  /*!
   * \brief specialization for constant fixed sized spans
   *
   * \tparam N: number of components of the span, which must be strictly
   * positive and known at compile time
   */
  template <std::size_t N>
  requires((N > 0) && (N != std::dynamic_extent))  //
      struct MaximumComponent<std::span<const real, N>> {
    /*!
     * \brief return the largest component of a span
     * \param[in] v: span holding `N` scalar values
     * \return the largest component of the given span
     * \note the comparisons are unrolled when `N` is lower than or equal to
     * four, `std::max_element` is used otherwise
     */
    [[nodiscard]] static constexpr real exe(const std::span<const real, N>& v)
        MGIS_NOEXCEPT
      requires((N != 0) && (N != std::dynamic_extent))
    {
      if constexpr (N == 1) {
        return v[0];
      } else if constexpr (N == 2) {
        return std::max(v[0], v[1]);
      } else if constexpr (N == 3) {
        return std::max(std::max(v[0], v[1]), v[2]);
      } else if constexpr (N == 4) {
        return std::max(std::max(std::max(v[0], v[1]), v[2]), v[3]);
      } else {
        return *(std::max_element(v.begin(), v.end()));
      }
    }
  };

  /*!
   * \brief specialization for fixed sized arrays of scalar values
   *
   * \tparam N: number of components of the array, which must be strictly
   * positive and known at compile time
   */
  template <std::size_t N>
  requires(N > 0) struct MaximumComponent<std::array<real, N>> {
    /*!
     * \brief return the largest component of an array
     * \param[in] v: array holding `N` scalar values
     * \return the largest component of the given array
     * \note the array is mapped to a constant span of the same size and the
     * call is forwarded to the specialization for spans
     */
    [[nodiscard]] static constexpr real exe(const std::array<real, N>& v)
        MGIS_NOEXCEPT {
      return MaximumComponent<std::span<const real, N>>::exe(
          std::span<const real, N>(v.data(), N));
    }
  };

  /*!
   * \brief customization point giving the smallest component of a value which
   * can be returned by an evaluator
   *
   * The primary template is only declared. Specializations are provided for the
   * result types supported by the library (`real`, fixed size constant spans
   * and arrays of scalar values) and new specializations can be added to
   * support other result types.
   *
   * \tparam T: type of the value whose smallest component is requested
   */
  template <typename T>
  struct MinimumComponent;

  /*!
   * \brief specialization for scalar values
   *
   * A scalar value has a single component: the smallest component is the value
   * itself.
   */
  template <>
  struct MinimumComponent<real> {
    /*!
     * \brief return the smallest component of a scalar value
     * \param[in] v: scalar value
     * \return the given value
     */
    [[nodiscard]] static constexpr real exe(const real& v) MGIS_NOEXCEPT {
      return v;
    }
  };

  /*!
   * \brief specialization for constant fixed sized spans
   *
   * \tparam N: number of components of the span, which must be strictly
   * positive and known at compile time
   */
  template <std::size_t N>
  requires((N > 0) && (N != std::dynamic_extent))  //
      struct MinimumComponent<std::span<const real, N>> {
    /*!
     * \brief return the smallest component of a span
     * \param[in] v: span holding `N` scalar values
     * \return the smallest component of the given span
     * \note the comparisons are unrolled when `N` is lower than or equal to
     * four, `std::min_element` is used otherwise
     */
    [[nodiscard]] static constexpr real exe(const std::span<const real, N>& v)
        MGIS_NOEXCEPT
      requires((N != 0) && (N != std::dynamic_extent))
    {
      if constexpr (N == 1) {
        return v[0];
      } else if constexpr (N == 2) {
        return std::min(v[0], v[1]);
      } else if constexpr (N == 3) {
        return std::min(std::min(v[0], v[1]), v[2]);
      } else if constexpr (N == 4) {
        return std::min(std::min(std::min(v[0], v[1]), v[2]), v[3]);
      } else {
        return *(std::min_element(v.begin(), v.end()));
      }
    }
  };

  /*!
   * \brief specialization for fixed sized arrays of scalar values
   *
   * \tparam N: number of components of the array, which must be strictly
   * positive and known at compile time
   */
  template <std::size_t N>
  requires(N > 0) struct MinimumComponent<std::array<real, N>> {
    /*!
     * \brief return the smallest component of an array
     * \param[in] v: array holding `N` scalar values
     * \return the smallest component of the given array
     * \note the array is mapped to a constant span of the same size and the
     * call is forwarded to the specialization for spans
     */
    [[nodiscard]] static constexpr real exe(const std::array<real, N>& v)
        MGIS_NOEXCEPT {
      return MinimumComponent<std::span<const real, N>>::exe(
          std::span<const real, N>(v.data(), N));
    }
  };

}  // namespace mgis::function::customization_points

namespace mgis::function::internals {

  /*!
   * \brief functor returning the absolute value of the values returned by an
   * evaluator, using the `customization_points::AbsoluteValue` customization
   * point
   */
  struct AbsoluteValueOperator {
    /*!
     * \brief call operator
     * \tparam ValueType: type of the value, which must be supported by the
     * `customization_points::AbsoluteValue` customization point
     * \param[in] v: value
     * \return the absolute value of each component of the given value
     */
    template <typename ValueType>
    MGIS_HOST_DEVICE [[nodiscard]] constexpr auto operator()(
        const ValueType& v) const {  //
      return customization_points::AbsoluteValue<std::decay_t<ValueType>>::exe(
          v);
    }
  };

  /*!
   * \brief functor returning the largest component of the values returned by an
   * evaluator, using the `customization_points::MaximumComponent` customization
   * point
   */
  struct MaximumComponentOperator {
    /*!
     * \brief call operator
     * \tparam ValueType: type of the value, which must be supported by the
     * `customization_points::MaximumComponent` customization point
     * \param[in] v: value
     * \return the largest component of the given value
     */
    template <typename ValueType>
    MGIS_HOST_DEVICE [[nodiscard]] constexpr auto operator()(
        const ValueType& v) const {  //
      return customization_points::MaximumComponent<
          std::decay_t<ValueType>>::exe(v);
    }
  };

  /*!
   * \brief functor returning the smallest component of the values returned by
   * an evaluator, using the `customization_points::MinimumComponent`
   * customization point
   */
  struct MinimumComponentOperator {
    /*!
     * \brief call operator
     * \tparam ValueType: type of the value, which must be supported by the
     * `customization_points::MinimumComponent` customization point
     * \param[in] v: value
     * \return the smallest component of the given value
     */
    template <typename ValueType>
    MGIS_HOST_DEVICE [[nodiscard]] constexpr auto operator()(
        const ValueType& v) const {  //
      return customization_points::MinimumComponent<
          std::decay_t<ValueType>>::exe(v);
    }
  };

#ifdef MGIS_HAVE_TFEL

  /*!
   * \brief functor returning the opposite of a value handled by `tfel::math`
   */
  struct NegateOperator {
    /*!
     * \brief call operator
     * \tparam OperandType: type of the operand, whose negation must be defined
     * in `tfel::math`
     * \param[in] a: operand
     * \return the opposite of the given operand, whose type is the result type
     * of the unary operation `tfel::math::OpNeg`
     * \note only available if the number of components of the result is known
     * at compile time
     */
    template <typename OperandType>
    MGIS_HOST_DEVICE [[nodiscard]] constexpr auto operator()(
        const OperandType& a) const
        -> tfel::math::UnaryOperationResult<OperandType,
                                            tfel::math::OpNeg>  //
    requires(compile_time_size<
                 tfel::math::UnaryOperationResult<OperandType,
                                                  tfel::math::OpNeg>> !=
             dynamic_extent) {  //
      return -a;
    }
  };

  /*!
   * \brief functor multiplying a value handled by `tfel::math` by a scalar
   * value
   */
  struct MultiplyByScalarOperator {
    /*!
     * \brief constructor
     * \param[in] s_: scalar value the operands are multiplied by
     */
    MGIS_HOST_DEVICE constexpr MultiplyByScalarOperator(const real s_)
        : s(s_) {}
    //
    /*!
     * \brief call operator
     * \tparam FirstOperandType: type of the operand, whose product by a scalar
     * value must be defined in `tfel::math`
     * \param[in] a: operand
     * \return the product of the operand by the stored scalar value, whose type
     * is the result type of the binary operation `tfel::math::OpMult`
     * \note only available if the number of components of the result is known
     * at compile time
     */
    template <typename FirstOperandType>
    MGIS_HOST_DEVICE [[nodiscard]] constexpr auto operator()(
        const FirstOperandType& a) const
        -> tfel::math::BinaryOperationResult<FirstOperandType,
                                             real,
                                             tfel::math::OpMult>  //
    requires(compile_time_size<
                 tfel::math::BinaryOperationResult<FirstOperandType,
                                                   real,
                                                   tfel::math::OpMult>> !=
             dynamic_extent) {
      return a * (this->s);
    }

   private:
    //! \brief scalar value the operands are multiplied by
    const real s;
  };

  /*!
   * \brief functor dividing a value handled by `tfel::math` by a scalar value
   */
  struct DivideByScalarOperator {
    /*!
     * \brief constructor
     * \param[in] s_: scalar value the operands are divided by
     */
    MGIS_HOST_DEVICE constexpr DivideByScalarOperator(const real s_) : s(s_) {}
    //
    /*!
     * \brief call operator
     * \tparam FirstOperandType: type of the operand, whose quotient by a scalar
     * value must be defined in `tfel::math`
     * \param[in] a: operand
     * \return the quotient of the operand by the stored scalar value, whose
     * type is the result type of the binary operation `tfel::math::OpMult`
     * \note the declared return type is the result type of the product by a
     * scalar value, which is the type of the quotient this operator returns
     * \note only available if the number of components of the result is known
     * at compile time
     */
    template <typename FirstOperandType>
    MGIS_HOST_DEVICE [[nodiscard]] constexpr auto operator()(
        const FirstOperandType& a) const
        -> tfel::math::BinaryOperationResult<FirstOperandType,
                                             real,
                                             tfel::math::OpMult>  //
    requires(compile_time_size<
                 tfel::math::BinaryOperationResult<FirstOperandType,
                                                   real,
                                                   tfel::math::OpMult>> !=
             dynamic_extent) {
      return a / (this->s);
    }

   private:
    //! \brief scalar value the operands are divided by
    const real s;
  };

#endif /* MGIS_HAVE_TFEL */

}  // end of namespace mgis::function::internals

namespace mgis::function {

  /*!
   * \brief modifier returning the absolute value of the values returned by an
   * evaluator, as in `e | absolute_value`
   *
   * The number of components of the values returned by the modifier is that of
   * the values computed by the underlying evaluator: the result is a scalar
   * value if the evaluator returns a scalar value, a `std::array` of scalar
   * values if it returns a fixed size span.
   */
  inline constexpr auto absolute_value =
      internals::unary_operation_modifier2(internals::AbsoluteValueOperator{});

  /*!
   * \brief modifier returning the largest component of the values returned by
   * an evaluator as a scalar value, as in `e | maximum_component`
   */
  inline constexpr auto maximum_component =
      internals::unary_operation_modifier2(
          internals::MaximumComponentOperator{});

  /*!
   * \brief modifier returning the smallest component of the values returned by
   * an evaluator as a scalar value, as in `e | minimum_component`
   */
  inline constexpr auto minimum_component =
      internals::unary_operation_modifier2(
          internals::MinimumComponentOperator{});

#ifdef MGIS_HAVE_TFEL

  /*!
   * \brief modifier returning the opposite of the values returned by an
   * evaluator, as in `e | negate`
   * \note this modifier requires the TFEL support: the number of components of
   * the opposite of the values returned by the evaluator must be known at
   * compile time
   */
  inline constexpr auto negate =
      internals::unary_operation_modifier2(internals::NegateOperator{});

  /*!
   * \brief create a modifier multiplying the values returned by an evaluator
   * by a scalar value, as in `e | multiply_by_scalar(2)`
   *
   * \param[in] s: scalar value the values are multiplied by
   * \return a modifier storing the scalar value
   * \note this modifier requires the TFEL support: the number of components of
   * the product must be known at compile time
   */
  constexpr auto multiply_by_scalar(const real s) {
    return internals::unary_operation_modifier<
        internals::MultiplyByScalarOperator>(
        internals::MultiplyByScalarOperator(s));
  }

  /*!
   * \brief create a modifier dividing the values returned by an evaluator by a
   * scalar value, as in `e | divide_by_scalar(3)`
   *
   * \param[in] s: scalar value the values are divided by
   * \return a modifier storing the scalar value
   * \note this modifier requires the TFEL support: the number of components of
   * the quotient must be known at compile time
   */
  constexpr auto divide_by_scalar(const real s) {
    return internals::unary_operation_modifier<
        internals::DivideByScalarOperator>(
        internals::DivideByScalarOperator(s));
  }  // end of divide_by_scalar

#endif /* MGIS_HAVE_TFEL */

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_UNARYOPERATION_HXX */
