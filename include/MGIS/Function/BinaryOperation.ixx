/*!
 * \file   MGIS/Function/BinaryOperation.ixx
 * \brief  Implementation of the `BinaryOperationModifier` and
 * `BinaryOperationModifier2` classes and of the helpers declared in
 * `MGIS/Function/BinaryOperation.hxx`.
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

#ifndef LIB_MGIS_FUNCTION_BINARYOPERATION_IXX
#define LIB_MGIS_FUNCTION_BINARYOPERATION_IXX

namespace mgis::function {

  // the callable object is copied in the object, the two evaluators are
  // handled by the `BinaryOperationEvaluatorBase` base class
  template <typename CallableType,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationModifierRequirements<CallableType,
                                               FirstEvaluatorType,
                                               SecondEvaluatorType>)  //
      MGIS_HOST_DEVICE constexpr BinaryOperationModifier<CallableType,
                                                         FirstEvaluatorType,
                                                         SecondEvaluatorType>::
          BinaryOperationModifier(const CallableType& c,
                                  const FirstEvaluatorType& e1,
                                  const SecondEvaluatorType& e2)
      : BinaryOperationEvaluatorBase<BinaryOperationModifier,
                                     FirstEvaluatorType,
                                     SecondEvaluatorType>(e1, e2),
        modifier(c) {}  // end of BinaryOperationModifier

  //! \brief calls the stored callable object with the given values
  template <typename CallableType,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationModifierRequirements<CallableType,
                                               FirstEvaluatorType,
                                               SecondEvaluatorType>)  //
      MGIS_HOST_DEVICE
      constexpr auto BinaryOperationModifier<CallableType,
                                             FirstEvaluatorType,
                                             SecondEvaluatorType>::
          apply(const evaluator_result<FirstEvaluatorType>& v1,
                const evaluator_result<SecondEvaluatorType>& v2) const {
    return this->modifier(v1, v2);
  }

  // the callable object is not stored in the object: a default
  // constructed object is used each time the values are computed, which is
  // allowed since `CallableType` is trivially default constructible
  template <typename CallableType,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationModifier2Requirements<CallableType,
                                                FirstEvaluatorType,
                                                SecondEvaluatorType>)  //
      MGIS_HOST_DEVICE
      constexpr auto BinaryOperationModifier2<CallableType,
                                              FirstEvaluatorType,
                                              SecondEvaluatorType>::
          apply(const evaluator_result<FirstEvaluatorType>& v1,
                const evaluator_result<SecondEvaluatorType>& v2) const {
    auto c = CallableType{};
    return c(v1, v2);
  }  // end of apply

  //! \brief forwards the call to the `getNumberOfComponents` member
  template <typename CallableType,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  MGIS_HOST_DEVICE constexpr mgis::size_type getNumberOfComponents(
      const BinaryOperationModifier<CallableType,
                                    FirstEvaluatorType,
                                    SecondEvaluatorType>& e) {
    return e.getNumberOfComponents();
  }  // end of getNumberOfComponents

  //! \brief forwards the call to the `getNumberOfComponents` member
  template <typename CallableType,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  MGIS_HOST_DEVICE constexpr mgis::size_type getNumberOfComponents(
      const BinaryOperationModifier2<CallableType,
                                     FirstEvaluatorType,
                                     SecondEvaluatorType>& e) {
    return e.getNumberOfComponents();
  }  // end of getNumberOfComponents

  namespace internals {

    //! \brief stores copies of the callable object and of the evaluator
    template <typename CallableType, EvaluatorConcept SecondEvaluatorType>
    requires(std::is_copy_constructible_v<CallableType>)  //
        MGIS_HOST_DEVICE
        constexpr BinaryOperatorCurrying<CallableType, SecondEvaluatorType>::
            BinaryOperatorCurrying(const CallableType& c,
                                   const SecondEvaluatorType& e)
        : modifier(c), e2(e) {}  // end of BinaryOperatorCurrying

    //! \brief builds a `BinaryOperationModifier` from the given evaluator
    template <typename CallableType, EvaluatorConcept SecondEvaluatorType>
    requires(std::is_copy_constructible_v<CallableType>)  //
        template <typename FirstEvaluatorType>
        MGIS_HOST_DEVICE constexpr auto BinaryOperatorCurrying<
            CallableType,
            SecondEvaluatorType>::operator()(FirstEvaluatorType&& e1) const
        requires((EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&(
            std::invocable<CallableType,
                           evaluator_result<std::decay_t<FirstEvaluatorType>>,
                           evaluator_result<SecondEvaluatorType>>)) {
      return BinaryOperationModifier<
          CallableType, std::decay_t<FirstEvaluatorType>, SecondEvaluatorType>{
          this->modifier, std::forward<FirstEvaluatorType>(e1), this->e2};
    }  // end of operator()

    //! \brief stores a copy of the evaluator of the second argument
    template <typename CallableType, EvaluatorConcept SecondEvaluatorType>
    requires(std::is_trivially_default_constructible_v<CallableType>)  //
        MGIS_HOST_DEVICE
        constexpr BinaryOperatorCurrying2<CallableType, SecondEvaluatorType>::
            BinaryOperatorCurrying2(const SecondEvaluatorType& e)
        : e2(e) {}

    //! \brief builds a `BinaryOperationModifier2` from the given evaluator
    template <typename CallableType, EvaluatorConcept SecondEvaluatorType>
    requires(std::is_trivially_default_constructible_v<CallableType>)  //
        template <typename FirstEvaluatorType>
        MGIS_HOST_DEVICE constexpr auto BinaryOperatorCurrying2<
            CallableType,
            SecondEvaluatorType>::operator()(FirstEvaluatorType&& e1) const
        requires(
            (EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&  //
            (std::invocable<CallableType,
                            evaluator_result<std::decay_t<FirstEvaluatorType>>,
                            evaluator_result<SecondEvaluatorType>>)) {
      return BinaryOperationModifier2<
          CallableType, std::decay_t<FirstEvaluatorType>, SecondEvaluatorType>(
          std::forward<FirstEvaluatorType>(e1), this->e2);
    }  // end of operator()

    //! \brief stores a copy of the callable object
    template <typename CallableType>
    MGIS_HOST_DEVICE constexpr binary_operation_modifier<
        CallableType>::binary_operation_modifier(const CallableType& c)
        : modifier(c) {}  // end of binary_operation_modifier

    //! \brief partially applies the operation to the given evaluator
    template <typename CallableType>
    template <typename SecondEvaluatorType>
    MGIS_HOST_DEVICE constexpr auto
    binary_operation_modifier<CallableType>::operator()(
        SecondEvaluatorType&& e2) const
        requires(EvaluatorConcept<std::decay_t<SecondEvaluatorType>>) {
      return BinaryOperatorCurrying<CallableType,
                                    std::decay_t<SecondEvaluatorType>>(
          this->modifier, std::forward<SecondEvaluatorType>(e2));
    }

    //! \brief builds a `BinaryOperationModifier` from the two evaluators
    template <typename CallableType>
    template <typename FirstEvaluatorType, typename SecondEvaluatorType>
    MGIS_HOST_DEVICE constexpr auto
    binary_operation_modifier<CallableType>::operator()(
        FirstEvaluatorType&& e1, SecondEvaluatorType&& e2) const
        requires((EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&   //
                 (EvaluatorConcept<std::decay_t<SecondEvaluatorType>>)&&  //
                 (std::invocable<
                     CallableType,
                     evaluator_result<std::decay_t<FirstEvaluatorType>>,
                     evaluator_result<std::decay_t<SecondEvaluatorType>>>)) {
      return BinaryOperationModifier<CallableType,
                                     std::decay_t<FirstEvaluatorType>,
                                     std::decay_t<SecondEvaluatorType>>(
          this->modifier, std::forward<FirstEvaluatorType>(e1),
          std::forward<SecondEvaluatorType>(e2));
    }  // end of operator()

    //! \brief partially applies the operation to the given evaluator
    template <typename CallableType>
    template <typename SecondEvaluatorType>
    MGIS_HOST_DEVICE constexpr auto
    binary_operation_modifier2_impl<CallableType>::operator()(
        SecondEvaluatorType&& e2) const
        requires(EvaluatorConcept<std::decay_t<SecondEvaluatorType>>) {
      return BinaryOperatorCurrying2<CallableType,
                                     std::decay_t<SecondEvaluatorType>>(
          std::forward<SecondEvaluatorType>(e2));
    }  // end of operator()

    //! \brief builds a `BinaryOperationModifier2` from the two evaluators
    template <typename CallableType>
    template <typename FirstEvaluatorType, typename SecondEvaluatorType>
    MGIS_HOST_DEVICE constexpr auto
    binary_operation_modifier2_impl<CallableType>::operator()(
        FirstEvaluatorType&& e1, SecondEvaluatorType&& e2) const
        requires((EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&   //
                 (EvaluatorConcept<std::decay_t<SecondEvaluatorType>>)&&  //
                 (std::invocable<
                     CallableType,
                     evaluator_result<std::decay_t<FirstEvaluatorType>>,
                     evaluator_result<std::decay_t<SecondEvaluatorType>>>)) {
      return BinaryOperationModifier2<CallableType,
                                      std::decay_t<FirstEvaluatorType>,
                                      std::decay_t<SecondEvaluatorType>>(
          std::forward<FirstEvaluatorType>(e1),
          std::forward<SecondEvaluatorType>(e2));
    }  // end of operator()

    //! \brief the argument is only used to deduce the type of the
    //!        binary operation
    template <typename CallableType>
    MGIS_HOST_DEVICE constexpr auto binary_operation_modifier2(CallableType) {
      return binary_operation_modifier2_impl<CallableType>{};
    }  // end of binary_operation_modifier2

  }  // namespace internals

  //! \brief creates a `binary_operation_modifier` from the callable object
  //!        and applies it to the given evaluator, which is the evaluator of
  //!        the second argument of the operation
  template <typename CallableType, typename SecondEvaluatorType>
  MGIS_HOST_DEVICE constexpr auto binary_operation(CallableType&& c,
                                                   SecondEvaluatorType&& e2)  //
      requires(EvaluatorConcept<std::decay_t<SecondEvaluatorType>>) {
    auto modifier =
        internals::binary_operation_modifier<std::decay_t<CallableType>>(
            std::forward<CallableType>(c));
    return modifier(std::forward<SecondEvaluatorType>(e2));
  }

  //! \brief applies the operation to the two given evaluators
  template <typename CallableType,
            typename FirstEvaluatorType,
            typename SecondEvaluatorType>
  constexpr auto binary_operation(CallableType&& c,
                                  FirstEvaluatorType&& e1,
                                  SecondEvaluatorType&& e2)             //
      requires((EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&   //
               (EvaluatorConcept<std::decay_t<SecondEvaluatorType>>)&&  //
               (std::invocable<
                   std::decay_t<CallableType>,
                   evaluator_result<std::decay_t<FirstEvaluatorType>>,
                   evaluator_result<std::decay_t<SecondEvaluatorType>>>)) {
    auto modifier =
        internals::binary_operation_modifier<std::decay_t<CallableType>>(
            std::forward<CallableType>(c));
    return modifier(std::forward<FirstEvaluatorType>(e1),
                    std::forward<SecondEvaluatorType>(e2));
  }  // end of binary_operation

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_BINARYOPERATION_IXX */
