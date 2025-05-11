/*!
 * \file   MGIS/Function/BinaryOperationEvaluator.ixx
 * \brief    
 * \author Thomas Helfer
 * \date   11/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATOR_IXX
#define LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATOR_IXX

namespace mgis::function {

  template <typename CallableType,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  BinaryOperationEvaluatorModifier<CallableType,
                                   FirstEvaluatorType,
                                   SecondEvaluatorType>::
      BinaryOperationEvaluatorModifier(const CallableType& c,
                                       const FirstEvaluatorType& e1,
                                       const SecondEvaluatorType& e2)
      : BinaryOperationEvaluatorBase<BinaryOperationEvaluatorModifier,
                                     FirstEvaluatorType,
                                     SecondEvaluatorType>(e1, e2),
        modifier(c) {}  // end of BinaryOperationEvaluatorModifier

  template <typename CallableType,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  auto BinaryOperationEvaluatorModifier<CallableType,
                                        FirstEvaluatorType,
                                        SecondEvaluatorType>::
      apply(const evaluator_result<FirstEvaluatorType>& v1,
            const evaluator_result<SecondEvaluatorType>& v2) const {
    return this->modifier(v1, v2);
  }

  template <typename CallableType,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  auto BinaryOperationEvaluatorModifier2<CallableType,
                                         FirstEvaluatorType,
                                         SecondEvaluatorType>::
      apply(const evaluator_result<FirstEvaluatorType>& v1,
            const evaluator_result<SecondEvaluatorType>& v2) const {
    auto c = CallableType{};
    return c(v1, v2);
  } // end of apply

  namespace internals {

    template <typename CallableType, EvaluatorConcept SecondEvaluatorType>
    BinaryOperatorCurrying<CallableType, SecondEvaluatorType>::
        BinaryOperatorCurrying(const CallableType& c,
                               const SecondEvaluatorType& e)
        : modifier(c), e2(e) {}  // end of BinaryOperatorCurrying

    template <typename CallableType, EvaluatorConcept SecondEvaluatorType>
    template <typename FirstEvaluatorType>
    auto BinaryOperatorCurrying<CallableType, SecondEvaluatorType>::operator()(
        FirstEvaluatorType&& e1) const
        requires((EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&(
            std::invocable<CallableType,
                           evaluator_result<std::decay_t<FirstEvaluatorType>>,
                           evaluator_result<SecondEvaluatorType>>)) {
      return BinaryOperationEvaluatorModifier(
          this->modifier, std::forward<FirstEvaluatorType>(e1), this->e2);
    }  // end of operator()

    template <typename CallableType, EvaluatorConcept SecondEvaluatorType>
    BinaryOperatorCurrying2<CallableType, SecondEvaluatorType>::
        BinaryOperatorCurrying2(const SecondEvaluatorType& e)
        : e2(e) {}

    template <typename CallableType, EvaluatorConcept SecondEvaluatorType>
    template <typename FirstEvaluatorType>
    auto BinaryOperatorCurrying2<CallableType, SecondEvaluatorType>::operator()(
        FirstEvaluatorType&& e1) const
        requires(
            (EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&  //
            (std::invocable<CallableType,
                            evaluator_result<std::decay_t<FirstEvaluatorType>>,
                            evaluator_result<SecondEvaluatorType>>)) {
      return BinaryOperationEvaluatorModifier2<
          CallableType, std::decay_t<FirstEvaluatorType>, SecondEvaluatorType>(
          std::forward<FirstEvaluatorType>(e1), this->e2);
    }  // end of operator()

    template <typename CallableType>
    binary_operation_modifier<CallableType>::binary_operation_modifier(
        CallableType&& c)
        : modifier(std::forward<CallableType>(c)) {
    }  // end of binary_operation_modifier

    template <typename CallableType>
    template <typename SecondEvaluatorType>
    auto binary_operation_modifier<CallableType>::operator()(
        SecondEvaluatorType&& e2) const
        requires(EvaluatorConcept<std::decay_t<SecondEvaluatorType>>) {
      return BinaryOperatorCurrying<CallableType,
                                    std::decay_t<SecondEvaluatorType>>(
          this->modifier, std::forward<SecondEvaluatorType>(e2));
    }

    template <typename CallableType>
    template <typename FirstEvaluatorType, typename SecondEvaluatorType>
    auto binary_operation_modifier<CallableType>::operator()(
        FirstEvaluatorType&& e1, SecondEvaluatorType&& e2) const
        requires((EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&   //
                 (EvaluatorConcept<std::decay_t<SecondEvaluatorType>>)&&  //
                 (std::invocable<
                     CallableType,
                     evaluator_result<std::decay_t<FirstEvaluatorType>>,
                     evaluator_result<std::decay_t<SecondEvaluatorType>>>)) {
      return BinaryOperationEvaluatorModifier<
          CallableType, std::decay_t<FirstEvaluatorType>,
          std::decay_t<SecondEvaluatorType>>(
          this->modifier, std::forward<FirstEvaluatorType>(e1),
          std::forward<SecondEvaluatorType>(e2));
      }  // end of operator()

      template <typename CallableType>
      template <typename SecondEvaluatorType>
      auto binary_operation_modifier2<CallableType>::operator()(
          SecondEvaluatorType&& e2) const
          requires(EvaluatorConcept<std::decay_t<SecondEvaluatorType>>) {
        return BinaryOperatorCurrying2<CallableType,
                                       std::decay_t<SecondEvaluatorType>>(
            std::forward<SecondEvaluatorType>(e2));
      }  // end of operator()

      template <typename CallableType>
      template <typename FirstEvaluatorType, typename SecondEvaluatorType>
      auto binary_operation_modifier2<CallableType>::operator()(
          FirstEvaluatorType&& e1, SecondEvaluatorType&& e2)
          const requires(
              (EvaluatorConcept<std::decay_t<FirstEvaluatorType>>)&&   //
              (EvaluatorConcept<std::decay_t<SecondEvaluatorType>>)&&  //
              (std::invocable<
                  CallableType,
                  evaluator_result<std::decay_t<FirstEvaluatorType>>,
                  evaluator_result<std::decay_t<SecondEvaluatorType>>>)) {
        return BinaryOperationEvaluatorModifier2<
            CallableType, std::decay_t<FirstEvaluatorType>,
            std::decay_t<SecondEvaluatorType>>(
            std::forward<FirstEvaluatorType>(e1),
            std::forward<SecondEvaluatorType>(e2));
      }  // end of operator()

    }  // namespace internals

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATOR_IXX */
