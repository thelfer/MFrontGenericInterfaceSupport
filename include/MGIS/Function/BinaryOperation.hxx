/*!
 * \file   MGIS/Function/BinaryOperation.hxx
 * \brief  This file declares the `BinaryOperationModifier` and
 * `BinaryOperationModifier2` classes, which build a new evaluator by
 * applying a binary operation to the values returned by two existing
 * evaluators, the `binary_operation` function used to build such
 * evaluators and the predefined modifiers `add`, `substract`, `multiply`,
 * `divide`, `mean_value` and `inner_product`.
 * \author Thomas Helfer
 * \date   11/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_BINARYOPERATION_HXX
#define LIB_MGIS_FUNCTION_BINARYOPERATION_HXX

#ifdef MGIS_HAVE_TFEL
#include "TFEL/Math/General/BasicOperations.hxx"
#include "TFEL/Math/General/ResultType.hxx"
#endif /* MGIS_HAVE_TFEL */

#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/CompileTimeSize.hxx"
#include "MGIS/Function/UnaryOperation.hxx"
#include "MGIS/Function/BinaryOperationEvaluatorBase.hxx"

namespace mgis::function {

  namespace internals {

    /*!
     * \brief helper base class giving the number of components of the
     * values returned by a binary operation applied to two evaluators
     *
     * The number of components is the compile-time size of the type
     * returned by the callable when it is called with the values returned
     * by the two evaluators. A static assertion is raised when this number
     * is only known at run time, as evaluators are required to know it at
     * compile time.
     *
     * \tparam CallableType: type of the binary operation
     * \tparam FirstEvaluatorType: type of the evaluator used to compute
     * the first argument of the operation
     * \tparam SecondEvaluatorType: type of the evaluator used to compute
     * the second argument of the operation
     */
    template <typename CallableType,
              EvaluatorConcept FirstEvaluatorType,
              EvaluatorConcept SecondEvaluatorType>
    struct BinaryOperationModifierBase {
      /*!
       * \brief return the number of components
       * \return number of components of the values returned by the binary
       * operation, which is known at compile time
       */
      MGIS_HOST_DEVICE [[nodiscard]] static constexpr size_type
      getNumberOfComponents() MGIS_NOEXCEPT {
        using result_type =
            std::invoke_result_t<CallableType,
                                 evaluator_result<FirstEvaluatorType>,
                                 evaluator_result<SecondEvaluatorType>>;
        return compile_time_size<result_type>;
      }
      static_assert(getNumberOfComponents() != dynamic_extent);
    };

  }  // end of namespace internals

  /*!
   * \brief requirements which must be satisfied by a binary operation
   * which is stored in the evaluator it applies to
   *
   * The callable object must be copy constructible and callable with the
   * values returned by the two evaluators it is applied to.
   *
   * \tparam CallableType: type of the binary operation
   * \tparam FirstEvaluatorType: type of the evaluator used to compute
   * the first argument of the operation
   * \tparam SecondEvaluatorType: type of the evaluator used to compute
   * the second argument of the operation
   */
  template <typename CallableType,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  inline constexpr auto BinaryOperationModifierRequirements =
      (std::is_copy_constructible_v<CallableType>)&&(
          std::invocable<CallableType,
                         evaluator_result<FirstEvaluatorType>,
                         evaluator_result<SecondEvaluatorType>>);

  /*!
   * \brief an evaluator applying a binary operation to the values
   * returned by two evaluators
   *
   * Such an evaluator is usually built with the `binary_operation`
   * function, which can either be applied to two evaluators at once or be
   * partially applied to the second evaluator only, in which case a
   * modifier is obtained and can be applied to the first evaluator using
   * `operator|`:
   *
   * \code{.cpp}
   * const auto e3 = binary_operation(c, e1, e2);
   * const auto e4 = e1 | binary_operation(c, e2);
   * \endcode
   *
   * The callable object is stored in the instance. The number of
   * components of the returned values is the compile-time size of the
   * type returned by the callable, which must thus be known at compile
   * time. The two evaluators must be defined on the same space and must
   * both be either element evaluators or quadrature evaluators, as stated
   * by `BinaryOperationEvaluatorBaseRequirement`. Only the call operators
   * which match the way both evaluators can be called are available.
   *
   * \tparam CallableType: type of the binary operation, which must be
   * copy constructible and callable with the values returned by the two
   * evaluators
   * \tparam FirstEvaluatorType: type of the evaluator used to compute
   * the first argument of the operation
   * \tparam SecondEvaluatorType: type of the evaluator used to compute
   * the second argument of the operation
   */
  template <typename CallableType,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationModifierRequirements<CallableType,
                                               FirstEvaluatorType,
                                               SecondEvaluatorType>)  //
      struct BinaryOperationModifier
      : BinaryOperationEvaluatorBase<
            BinaryOperationModifier<CallableType,
                                    FirstEvaluatorType,
                                    SecondEvaluatorType>,
            FirstEvaluatorType,
            SecondEvaluatorType>,
        internals::BinaryOperationModifierBase<CallableType,
                                               FirstEvaluatorType,
                                               SecondEvaluatorType> {
    /*!
     * \brief constructor
     * \param[in] c: binary operation, a copy of which is stored in the object
     * \param[in] e1: evaluator used to compute the first argument of the
     * operation
     * \param[in] e2: evaluator used to compute the second argument of the
     * operation
     */
    MGIS_HOST_DEVICE constexpr BinaryOperationModifier(
        const CallableType&,
        const FirstEvaluatorType&,
        const SecondEvaluatorType&);
    //
    using internals::BinaryOperationModifierBase<
        CallableType,
        FirstEvaluatorType,
        SecondEvaluatorType>::getNumberOfComponents;
    /*!
     * \brief apply the modifier
     * \param[in] v1: value computed by the evaluator of the first argument of
     * the operation at a given integration point
     * \param[in] v2: value computed by the evaluator of the second argument of
     * the operation at the same integration point
     * \return values computed by calling the stored callable object with
     * the values `v1` and `v2` computed by the evaluators of the first and
     * second argument of the operation at a given integration point
     * \note this function is called by the call operators inherited from
     * `BinaryOperationEvaluatorBase`
     */
    MGIS_HOST_DEVICE [[nodiscard]] constexpr auto apply(
        const evaluator_result<FirstEvaluatorType>&,
        const evaluator_result<SecondEvaluatorType>&) const;

   private:
    //! \brief binary operation applied to the values of the evaluators
    CallableType modifier;
  };

  /*!
   * \brief return the number of components of the values returned by a
   * binary operation
   * \tparam CallableType: type of the binary operation
   * \tparam FirstEvaluatorType: type of the evaluator used to compute
   * the first argument of the operation
   * \tparam SecondEvaluatorType: type of the evaluator used to compute
   * the second argument of the operation
   * \param[in] e: binary operation
   * \return number of components of the values returned by the operation
   */
  template <typename CallableType,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  MGIS_HOST_DEVICE [[nodiscard]] constexpr mgis::size_type
  getNumberOfComponents(const BinaryOperationModifier<CallableType,
                                                      FirstEvaluatorType,
                                                      SecondEvaluatorType>&);

  /*!
   * \brief requirements which must be satisfied by a binary operation
   * which is not stored in the evaluator it applies to
   *
   * The callable object must be trivially default constructible and
   * callable with the values returned by the two evaluators it is applied to.
   *
   * \tparam CallableType: type of the binary operation
   * \tparam FirstEvaluatorType: type of the evaluator used to compute
   * the first argument of the operation
   * \tparam SecondEvaluatorType: type of the evaluator used to compute
   * the second argument of the operation
   */
  template <typename CallableType,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  inline constexpr auto BinaryOperationModifier2Requirements =
      (std::is_trivially_default_constructible_v<CallableType>)&&(
          std::invocable<CallableType,
                         evaluator_result<FirstEvaluatorType>,
                         evaluator_result<SecondEvaluatorType>>);

  /*!
   * \brief an evaluator applying a binary operation to the values
   * returned by two evaluators, without storing this operation
   *
   * This class behaves as `BinaryOperationModifier` except that the
   * callable object is not stored in the instance: it must be trivially
   * default constructible and a default constructed object is created
   * each time the values are computed. Such a modifier is thus as light
   * as the two underlying evaluators, which allows `constexpr` modifiers
   * such as `add` or `inner_product` to be defined:
   *
   * \code{.cpp}
   * const auto e3 = add(e1, e2);
   * const auto e4 = e1 | add(e2);
   * \endcode
   *
   * The constructor inherited from `BinaryOperationEvaluatorBase` only
   * takes the two evaluators as arguments.
   *
   * \tparam CallableType: type of the binary operation, which must be
   * trivially default constructible and callable with the values returned
   * by the two evaluators
   * \tparam FirstEvaluatorType: type of the evaluator used to compute
   * the first argument of the operation
   * \tparam SecondEvaluatorType: type of the evaluator used to compute
   * the second argument of the operation
   */
  template <typename CallableType,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationModifier2Requirements<CallableType,
                                                FirstEvaluatorType,
                                                SecondEvaluatorType>)  //
      struct BinaryOperationModifier2
      : BinaryOperationEvaluatorBase<
            BinaryOperationModifier2<CallableType,
                                     FirstEvaluatorType,
                                     SecondEvaluatorType>,
            FirstEvaluatorType,
            SecondEvaluatorType>,
        internals::BinaryOperationModifierBase<CallableType,
                                               FirstEvaluatorType,
                                               SecondEvaluatorType> {
    // inheriting constructor
    using BinaryOperationEvaluatorBase<
        BinaryOperationModifier2,
        FirstEvaluatorType,
        SecondEvaluatorType>::BinaryOperationEvaluatorBase;
    //
    using internals::BinaryOperationModifierBase<
        CallableType,
        FirstEvaluatorType,
        SecondEvaluatorType>::getNumberOfComponents;
    /*!
     * \brief apply the modifier
     * \param[in] v1: value computed by the evaluator of the first argument of
     * the operation at a given integration point
     * \param[in] v2: value computed by the evaluator of the second argument of
     * the operation at the same integration point
     * \return values computed by calling a default constructed object of
     * type `CallableType` with the values `v1` and `v2` computed by the
     * evaluators of the first and second argument of the operation at a
     * given integration point
     * \note this function is called by the call operators inherited from
     * `BinaryOperationEvaluatorBase`
     */
    MGIS_HOST_DEVICE constexpr auto apply(
        const evaluator_result<FirstEvaluatorType>&,
        const evaluator_result<SecondEvaluatorType>&) const;
  };

  /*!
   * \brief return the number of components of the values returned by a
   * binary operation which does not store its callable object
   * \tparam CallableType: type of the binary operation
   * \tparam FirstEvaluatorType: type of the evaluator used to compute
   * the first argument of the operation
   * \tparam SecondEvaluatorType: type of the evaluator used to compute
   * the second argument of the operation
   * \param[in] e: binary operation
   * \return number of components of the values returned by the operation
   */
  template <typename CallableType,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  MGIS_HOST_DEVICE [[nodiscard]] constexpr mgis::size_type
  getNumberOfComponents(const BinaryOperationModifier2<CallableType,
                                                       FirstEvaluatorType,
                                                       SecondEvaluatorType>&);

  namespace internals {

    /*!
     * \brief partial application of a binary operation to the evaluator of
     * its second argument
     *
     * This helper matches the `EvaluatorModifierConcept` concept: applying
     * it to the evaluator of the first argument, with `operator|` for
     * instance, gives the evaluator of the result of the binary operation.
     * Both the callable object and the evaluator of the second argument are
     * stored in the object.
     *
     * \tparam CallableType: type of the binary operation, which must be
     * copy constructible
     * \tparam SecondEvaluatorType: type of the evaluator used to compute
     * the second argument of the operation
     */
    template <typename CallableType, EvaluatorConcept SecondEvaluatorType>
    requires(std::is_copy_constructible_v<
             CallableType>) struct BinaryOperatorCurrying {
      //! \brief this alias allows to match the Evaluator Modifier concept
      using Tag = ::mgis::function::EvaluatorModifierTag;
      //
      /*!
       * \brief constructor
       * \param[in] c: binary operation, a copy of which is stored in the object
       * \param[in] e: evaluator of the second argument of the operation, a
       * copy of which is stored in the object
       */
      MGIS_HOST_DEVICE constexpr BinaryOperatorCurrying(
          const CallableType&, const SecondEvaluatorType&);
      //
      /*!
       * \brief apply the partially applied binary operation to the
       * evaluator of its first argument
       * \tparam FirstEvaluatorType: type of the evaluator, which may be a
       * qualified type
       * \param[in] e1: evaluator used to compute the first argument of the
       * operation
       * \return a `BinaryOperationModifier` object, i.e. an evaluator
       * returning the values computed by the binary operation called with
       * the values of the evaluator `e1` used to compute the first argument
       */
      template <typename FirstEvaluatorType>
      MGIS_HOST_DEVICE constexpr auto
      operator()(FirstEvaluatorType&&) const requires(
          (EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&  //
          (std::invocable<CallableType,
                          evaluator_result<std::decay_t<FirstEvaluatorType>>,
                          evaluator_result<SecondEvaluatorType>>));

     private:
      //! \brief binary operation applied to the values of the evaluators
      CallableType modifier;
      //! \brief evaluator of the second argument of the binary operation
      SecondEvaluatorType e2;
    };

    /*!
     * \brief partial application of a binary operation to the evaluator of
     * its second argument, without storing the operation
     *
     * This helper behaves as `BinaryOperatorCurrying` except that the
     * callable object is not stored in the object: it must be trivially
     * default constructible. Only the evaluator of the second argument is
     * stored, which allows the predefined modifiers of this file to be
     * defined as `constexpr` global variables.
     *
     * \tparam CallableType: type of the binary operation, which must be
     * trivially default constructible
     * \tparam SecondEvaluatorType: type of the evaluator used to compute
     * the second argument of the operation
     */
    template <typename CallableType, EvaluatorConcept SecondEvaluatorType>
    requires(std::is_trivially_default_constructible_v<
             CallableType>) struct BinaryOperatorCurrying2 {
      //! \brief this alias allows to match the Evaluator Modifier concept
      using Tag = ::mgis::function::EvaluatorModifierTag;
      //
      /*!
       * \brief constructor
       * \param[in] e: evaluator of the second argument of the operation, a
       * copy of which is stored in the object
       */
      MGIS_HOST_DEVICE constexpr BinaryOperatorCurrying2(
          const SecondEvaluatorType&);
      //
      /*!
       * \brief apply the partially applied binary operation to the
       * evaluator of its first argument
       * \tparam FirstEvaluatorType: type of the evaluator, which may be a
       * qualified type
       * \param[in] e1: evaluator used to compute the first argument of the
       * operation
       * \return a `BinaryOperationModifier2` object, i.e. an evaluator
       * returning the values computed by the binary operation called with
       * the values of the evaluator `e1` used to compute the first argument
       */
      template <typename FirstEvaluatorType>
      MGIS_HOST_DEVICE constexpr auto
      operator()(FirstEvaluatorType&&) const requires(
          (EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&  //
          (std::invocable<CallableType,
                          evaluator_result<std::decay_t<FirstEvaluatorType>>,
                          evaluator_result<SecondEvaluatorType>>));

     private:
      //! \brief evaluator of the second argument of the binary operation
      SecondEvaluatorType e2;
    };

    /*!
     * \brief helper applying a stored binary operation to one or two evaluators
     *
     * The binary operation is stored in the helper. Calling the helper with
     * the evaluator of the second argument only gives a
     * `BinaryOperatorCurrying` object, i.e. a modifier which is applied to
     * the evaluator of the first argument. Calling the helper with the two
     * evaluators gives the evaluator of the result of the binary operation.
     *
     * \tparam CallableType: type of the binary operation, which must be
     * copy constructible
     * \note this helper does not match the `EvaluatorModifierConcept`
     * concept itself: its single argument call operator returns such a modifier
     */
    template <typename CallableType>
    struct binary_operation_modifier {
      //
      /*!
       * \brief constructor
       * \param[in] c: binary operation, a copy of it is stored in the helper
       */
      MGIS_HOST_DEVICE constexpr binary_operation_modifier(const CallableType&);
      //
      /*!
       * \brief partially apply the binary operation to the evaluator of
       * its second argument
       * \tparam SecondEvaluatorType: type of the evaluator, which may be
       * a qualified type
       * \param[in] e2: evaluator used to compute the second argument
       * \return a `BinaryOperatorCurrying` object, which can be applied to
       * the evaluator of the first argument using `operator|`
       */
      template <typename SecondEvaluatorType>
      MGIS_HOST_DEVICE constexpr auto operator()(SecondEvaluatorType&&) const
          requires(EvaluatorConcept<std::decay_t<SecondEvaluatorType>>);
      //
      /*!
       * \brief build the evaluator of the result of the binary operation
       * \tparam FirstEvaluatorType: type of the evaluator of the first
       * argument, which may be a qualified type
       * \tparam SecondEvaluatorType: type of the evaluator of the second
       * argument, which may be a qualified type
       * \param[in] e1: evaluator used to compute the first argument
       * \param[in] e2: evaluator used to compute the second argument
       * \return a `BinaryOperationModifier` object, i.e. an evaluator
       * returning the values computed by the binary operation
       */
      template <typename FirstEvaluatorType, typename SecondEvaluatorType>
      MGIS_HOST_DEVICE constexpr auto operator()(FirstEvaluatorType&&,
                                                 SecondEvaluatorType&&) const
          requires((EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&   //
                   (EvaluatorConcept<std::decay_t<SecondEvaluatorType>>)&&  //
                   (std::invocable<
                       CallableType,
                       evaluator_result<std::decay_t<FirstEvaluatorType>>,
                       evaluator_result<std::decay_t<SecondEvaluatorType>>>));

     private:
      //! \brief binary operation applied to the values of the evaluators
      CallableType modifier;
    };

    /*!
     * \brief helper applying a binary operation which is not stored in the
     * helper to one or two evaluators
     *
     * This helper behaves as `binary_operation_modifier` except that it is a
     * stateless type: the binary operation is created when the values are
     * computed, which allows the predefined modifiers of this file to be
     * defined as `constexpr` global variables.
     *
     * \tparam CallableType: type of the binary operation, which must be
     * trivially default constructible
     * \note this helper does not match the `EvaluatorModifierConcept`
     * concept itself: its single argument call operator returns such a modifier
     */
    template <typename CallableType>
    struct binary_operation_modifier2_impl {
      //
      /*!
       * \brief partially apply the binary operation to the evaluator of
       * its second argument
       * \tparam SecondEvaluatorType: type of the evaluator, which may be
       * a qualified type
       * \param[in] e2: evaluator used to compute the second argument
       * \return a `BinaryOperatorCurrying2` object, which can be applied
       * to the evaluator of the first argument using `operator|`
       */
      template <typename SecondEvaluatorType>
      MGIS_HOST_DEVICE constexpr auto operator()(SecondEvaluatorType&&) const
          requires(EvaluatorConcept<std::decay_t<SecondEvaluatorType>>);
      //
      /*!
       * \brief build the evaluator of the result of the binary operation
       * \tparam FirstEvaluatorType: type of the evaluator of the first
       * argument, which may be a qualified type
       * \tparam SecondEvaluatorType: type of the evaluator of the second
       * argument, which may be a qualified type
       * \param[in] e1: evaluator used to compute the first argument
       * \param[in] e2: evaluator used to compute the second argument
       * \return a `BinaryOperationModifier2` object, i.e. an evaluator
       * returning the values computed by the binary operation
       */
      template <typename FirstEvaluatorType, typename SecondEvaluatorType>
      MGIS_HOST_DEVICE constexpr auto operator()(FirstEvaluatorType&&,
                                                 SecondEvaluatorType&&) const
          requires((EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&   //
                   (EvaluatorConcept<std::decay_t<SecondEvaluatorType>>)&&  //
                   (std::invocable<
                       CallableType,
                       evaluator_result<std::decay_t<FirstEvaluatorType>>,
                       evaluator_result<std::decay_t<SecondEvaluatorType>>>));
    };

    /*!
     * \brief create a helper from a callable object which is only used to
     * deduce the type of the binary operation
     *
     * \tparam CallableType: type of the binary operation, which must be
     * trivially default constructible
     * \param[in] c: callable object which is discarded, as the operation is
     * recreated when the values are computed
     * \return a stateless helper of type
     * `binary_operation_modifier2_impl<CallableType>`
     */
    template <typename CallableType>
    MGIS_HOST_DEVICE constexpr auto binary_operation_modifier2(CallableType);

  }  // namespace internals

  /*!
   * \brief partially apply a binary operation to the evaluator of its
   * second argument
   *
   * The returned modifier can be applied to the evaluator of the first
   * argument using `operator|`, as in `e1 | binary_operation(c, e2)`.
   *
   * \tparam CallableType: type of the binary operation, which may be a
   * qualified type
   * \tparam SecondEvaluatorType: type of the evaluator of the second
   * argument, which may be a qualified type
   * \param[in] c: binary operation, a copy of it is stored in the modifier
   * \param[in] e2: evaluator used to compute the second argument
   * \return a `internals::BinaryOperatorCurrying` object
   */
  template <typename CallableType, typename SecondEvaluatorType>
  MGIS_HOST_DEVICE constexpr auto binary_operation(CallableType&&,
                                                   SecondEvaluatorType&&)  //
      requires(EvaluatorConcept<std::decay_t<SecondEvaluatorType>>);

  /*!
   * \brief apply a binary operation to the values returned by two evaluators
   *
   * \tparam CallableType: type of the binary operation, which may be a
   * qualified type
   * \tparam FirstEvaluatorType: type of the evaluator of the first
   * argument, which may be a qualified type
   * \tparam SecondEvaluatorType: type of the evaluator of the second
   * argument, which may be a qualified type
   * \param[in] c: binary operation, which must be callable with the values
   * returned by the two evaluators
   * \param[in] e1: evaluator used to compute the first argument
   * \param[in] e2: evaluator used to compute the second argument
   * \return a `BinaryOperationModifier` object, i.e. an evaluator returning
   * the values computed by the binary operation
   */
  template <typename CallableType,
            typename FirstEvaluatorType,
            typename SecondEvaluatorType>
  MGIS_HOST_DEVICE constexpr auto binary_operation(CallableType&&,
                                                   FirstEvaluatorType&&,
                                                   SecondEvaluatorType&&)  //
      requires((EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&      //
               (EvaluatorConcept<std::decay_t<SecondEvaluatorType>>)&&     //
               (std::invocable<
                   std::decay_t<CallableType>,
                   evaluator_result<std::decay_t<FirstEvaluatorType>>,
                   evaluator_result<std::decay_t<SecondEvaluatorType>>>));

}  // end of namespace mgis::function

#include "MGIS/Function/BinaryOperation.ixx"

#ifdef MGIS_HAVE_TFEL

namespace mgis::function::internals {

  /*!
   * \brief functor computing the sum of two values handled by `tfel::math`
   */
  struct AdditionOperator {
    //
    /*!
     * \brief call operator
     * \tparam FirstOperandType: type of the first operand, whose sum with
     * the second operand must be defined in `tfel::math`
     * \tparam SecondOperandType: type of the second operand
     * \param[in] a: first operand
     * \param[in] b: second operand
     * \return the sum of the two operands, whose type is the result type
     * of the binary operation `tfel::math::OpPlus`
     * \note only available if the number of components of the result is
     * known at compile time
     */
    template <typename FirstOperandType, typename SecondOperandType>
    MGIS_HOST_DEVICE [[nodiscard]] constexpr auto operator()(
        const FirstOperandType& a, const SecondOperandType& b) const
        -> tfel::math::BinaryOperationResult<FirstOperandType,
                                             SecondOperandType,
                                             tfel::math::OpPlus>  //
    requires(compile_time_size<
                 tfel::math::BinaryOperationResult<FirstOperandType,
                                                   SecondOperandType,
                                                   tfel::math::OpPlus>> !=
             dynamic_extent) {
      return a + b;
    }
  };

  /*!
   * \brief functor computing the difference of two values handled by
   * `tfel::math`
   */
  struct SubstractionOperator {
    //
    /*!
     * \brief call operator
     * \tparam FirstOperandType: type of the first operand, whose
     * difference with the second operand must be defined in `tfel::math`
     * \tparam SecondOperandType: type of the second operand
     * \param[in] a: first operand, which is the first term of the difference
     * \param[in] b: second operand, which is subtracted from the first one
     * \return the difference of the two operands, whose type is the result
     * type of the binary operation `tfel::math::OpMinus`
     * \note only available if the number of components of the result is
     * known at compile time
     * \note the call operator is not declared `const`, which is possible
     * since the functor has no data member
     */
    template <typename FirstOperandType, typename SecondOperandType>
    MGIS_HOST_DEVICE constexpr auto operator()(const FirstOperandType& a,
                                               const SecondOperandType& b)
        -> tfel::math::BinaryOperationResult<FirstOperandType,
                                             SecondOperandType,
                                             tfel::math::OpMinus>  //
    requires(compile_time_size<
                 tfel::math::BinaryOperationResult<FirstOperandType,
                                                   SecondOperandType,
                                                   tfel::math::OpMinus>> !=
             dynamic_extent) {
      return a - b;
    }
  };

  /*!
   * \brief functor computing the product of two values handled by
   * `tfel::math`
   */
  struct MultiplyOperator {
    //
    /*!
     * \brief call operator
     * \tparam FirstOperandType: type of the first operand, whose product
     * with the second operand must be defined in `tfel::math`
     * \tparam SecondOperandType: type of the second operand
     * \param[in] a: first operand
     * \param[in] b: second operand
     * \return the product of the two operands, whose type is the result
     * type of the binary operation `tfel::math::OpMult`
     * \note only available if the number of components of the result is
     * known at compile time
     */
    template <typename FirstOperandType, typename SecondOperandType>
    MGIS_HOST_DEVICE [[nodiscard]] constexpr auto operator()(
        const FirstOperandType& a, const SecondOperandType& b) const
        -> tfel::math::BinaryOperationResult<FirstOperandType,
                                             SecondOperandType,
                                             tfel::math::OpMult>  //
    requires(compile_time_size<
                 tfel::math::BinaryOperationResult<FirstOperandType,
                                                   SecondOperandType,
                                                   tfel::math::OpMult>> !=
             dynamic_extent) {
      return a * b;
    }
  };

  /*!
   * \brief functor computing the quotient of two values handled by
   * `tfel::math`
   */
  struct DivideOperator {
    //
    /*!
     * \brief call operator
     * \tparam FirstOperandType: type of the first operand, whose quotient
     * by the second operand must be defined in `tfel::math`
     * \tparam SecondOperandType: type of the second operand
     * \param[in] a: first operand, which is the dividend
     * \param[in] b: second operand, which is the divisor
     * \return the quotient of the two operands, whose type is the result
     * type of the binary operation `tfel::math::OpDiv`
     * \note only available if the number of components of the result is
     * known at compile time
     */
    template <typename FirstOperandType, typename SecondOperandType>
    MGIS_HOST_DEVICE [[nodiscard]] constexpr auto operator()(
        const FirstOperandType& a, const SecondOperandType& b) const
        -> tfel::math::BinaryOperationResult<FirstOperandType,
                                             SecondOperandType,
                                             tfel::math::OpDiv>  //
    requires(compile_time_size<
                 tfel::math::BinaryOperationResult<FirstOperandType,
                                                   SecondOperandType,
                                                   tfel::math::OpDiv>> !=
             dynamic_extent) {
      return a / b;
    }
  };

  /*!
   * \brief functor computing the mean value of two values handled by
   * `tfel::math`
   */
  struct MeanValueOperator {
    //
    /*!
     * \brief call operator
     * \tparam FirstOperandType: type of the first operand, whose mean
     * value with the second operand must be defined in `tfel::math`
     * \tparam SecondOperandType: type of the second operand
     * \param[in] a: first operand
     * \param[in] b: second operand
     * \return the expression `(a + b) / 2`, whose declared type is the result
     * type of the quotient by a scalar value of the result type of the binary
     * operation `tfel::math::OpMinus`
     * \note only available if the number of components of the result is
     * known at compile time
     */
    template <typename FirstOperandType, typename SecondOperandType>
    MGIS_HOST_DEVICE [[nodiscard]] constexpr auto operator()(
        const FirstOperandType& a, const SecondOperandType& b) const
        -> tfel::math::BinaryOperationResult<
            tfel::math::BinaryOperationResult<FirstOperandType,
                                              SecondOperandType,
                                              tfel::math::OpMinus>,
            real,
            tfel::math::OpDiv>  //
    requires(compile_time_size<tfel::math::BinaryOperationResult<
                 tfel::math::BinaryOperationResult<FirstOperandType,
                                                   SecondOperandType,
                                                   tfel::math::OpMinus>,
                 real,
                 tfel::math::OpDiv>> != dynamic_extent) {
      return (a + b) / 2;
    }
  };

  /*!
   * \brief functor computing the inner product of two values handled by
   * `tfel::math`
   */
  struct InnerProductOperator {
    //
    /*!
     * \brief call operator
     * \tparam FirstOperandType: type of the first operand, whose inner
     * product with the second operand must be defined in `tfel::math`
     * \tparam SecondOperandType: type of the second operand
     * \param[in] a: first operand
     * \param[in] b: second operand
     * \return the inner product of the two operands, whose type is the
     * result type of the binary operation `tfel::math::OpDotProduct`
     * \note only available if the number of components of the result is
     * known at compile time
     */
    template <typename FirstOperandType, typename SecondOperandType>
    MGIS_HOST_DEVICE [[nodiscard]] constexpr auto operator()(
        const FirstOperandType& a, const SecondOperandType& b) const
        -> tfel::math::BinaryOperationResult<FirstOperandType,
                                             SecondOperandType,
                                             tfel::math::OpDotProduct>  //
    requires(compile_time_size<
                 tfel::math::BinaryOperationResult<FirstOperandType,
                                                   SecondOperandType,
                                                   tfel::math::OpDotProduct>> !=
             dynamic_extent) {
      return a | b;
    }
  };

}  // end of namespace mgis::function::internals

namespace mgis::function {

  /*!
   * \brief predefined binary operation adding the values returned by two
   * evaluators: `add(e2)` gives a modifier which is applied to the evaluator
   * of the first argument, as in `e1 | add(e2)`, whereas `add(e1, e2)` gives
   * directly the evaluator of the sum
   * \note this operation requires the TFEL support
   */
  inline constexpr auto add =
      internals::binary_operation_modifier2(internals::AdditionOperator{});

  /*!
   * \brief predefined binary operation subtracting the values returned by an
   * evaluator from those returned by another one, as in `e1 | substract(e2)`
   * or `substract(e1, e2)`
   * \note this operation requires the TFEL support
   */
  inline constexpr auto substract =
      internals::binary_operation_modifier2(internals::SubstractionOperator{});

  /*!
   * \brief predefined binary operation multiplying the values returned by two
   * evaluators, as in `e1 | multiply(e2)` or `multiply(e1, e2)`
   * \note this operation requires the TFEL support
   */
  inline constexpr auto multiply =
      internals::binary_operation_modifier2(internals::MultiplyOperator{});

  /*!
   * \brief predefined binary operation dividing the values returned by an
   * evaluator by those returned by another one, as in `e1 | divide(e2)` or
   * `divide(e1, e2)`
   * \note this operation requires the TFEL support
   */
  inline constexpr auto divide =
      internals::binary_operation_modifier2(internals::DivideOperator{});

  /*!
   * \brief predefined binary operation computing the mean value of the values
   * returned by two evaluators, as in `e1 | mean_value(e2)` or
   * `mean_value(e1, e2)`
   * \note this operation requires the TFEL support
   */
  inline constexpr auto mean_value =
      internals::binary_operation_modifier2(internals::MeanValueOperator{});

  /*!
   * \brief predefined binary operation computing the inner product of the
   * values returned by two evaluators, which gives scalar values, as in
   * `e1 | inner_product(e2)` or `inner_product(e1, e2)`
   * \note this operation requires the TFEL support
   */
  inline constexpr auto inner_product =
      internals::binary_operation_modifier2(internals::InnerProductOperator{});

}  // end of namespace mgis::function

#endif /* MGIS_HAVE_TFEL */

#endif /* LIB_MGIS_FUNCTIONBINARYOPERATION_HXX */
