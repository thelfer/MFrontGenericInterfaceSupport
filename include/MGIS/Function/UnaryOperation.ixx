/*!
 * \file   MGIS/Function/UnaryOperation.ixx
 * \brief
 * \author Thomas Helfer
 * \date   09/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_UNARYOPERATION_IXX
#define LIB_MGIS_FUNCTION_UNARYOPERATION_IXX

namespace mgis::function {

  template <typename CallableType, EvaluatorConcept EvaluatorType>
  requires((std::is_copy_constructible_v<CallableType>)&&(
      std::invocable<CallableType,
                     evaluator_result<EvaluatorType>>))  //
    constexpr      UnaryOperation<CallableType, EvaluatorType>::UnaryOperation(
          const CallableType& c, const EvaluatorType& e)
      : EvaluatorModifierBase<UnaryOperation<CallableType, EvaluatorType>,
                              EvaluatorType>(e),
        modifier(c) {}  // end of UnaryOperation

  template <typename CallableType, EvaluatorConcept EvaluatorType>
  requires((std::is_copy_constructible_v<CallableType>)&&(
      std::invocable<CallableType,
                     evaluator_result<EvaluatorType>>))  //
    constexpr      auto UnaryOperation<CallableType, EvaluatorType>::apply(
          const evaluator_result<EvaluatorType>& values) const {
    return this->modifier(values);
  }  // end of apply

  template <typename CallableType, EvaluatorConcept EvaluatorType>
  requires((std::is_trivially_default_constructible_v<CallableType>)&&(
      std::invocable<CallableType,
                     evaluator_result<EvaluatorType>>))  //
    constexpr      auto UnaryOperation2<CallableType, EvaluatorType>::apply(
          const evaluator_result<EvaluatorType>& values) const {
    auto c = CallableType{};
    return c(values);
  }  // end of apply

  namespace internals {

    template <typename CallableType>
    constexpr    unary_operation_modifier<CallableType>::unary_operation_modifier(
        const CallableType& c)
        : modifier(c) {}

    template <typename CallableType>
    template <typename EvaluatorType>
    constexpr    auto unary_operation_modifier<CallableType>::operator()(EvaluatorType&& e)
        const requires((EvaluatorConcept<std::decay_t<EvaluatorType>>)&&(
            std::invocable<CallableType,
                           evaluator_result<std::decay_t<EvaluatorType>>>)) {
      return UnaryOperation<CallableType, std::decay_t<EvaluatorType>>(
          this->modifier, std::forward<EvaluatorType>(e));
    }  // end of operator()

    template <typename CallableType>
    template <typename EvaluatorType>
    constexpr    auto unary_operation_modifier2_impl<CallableType>::operator()(
        EvaluatorType&& e) const
        requires((EvaluatorConcept<std::decay_t<EvaluatorType>>)&&(
            std::invocable<CallableType,
                           evaluator_result<std::decay_t<EvaluatorType>>>)) {
      return UnaryOperation2<CallableType, std::decay_t<EvaluatorType>>(
          std::forward<EvaluatorType>(e));
    }  // end of operator()

    template <typename CallableType>
    constexpr auto unary_operation_modifier2(CallableType) {
      return unary_operation_modifier2_impl<CallableType>{};
    }  // end of unary_operation_modifier2

  }  // namespace internals

  template <typename CallableType>
  constexpr  auto unary_operation(CallableType&& c) {
    return internals::unary_operation_modifier<std::decay_t<CallableType>>(
        std::forward<CallableType>(c));
  }  // end of unary_operation

  template <typename CallableType>
  constexpr  auto transform(CallableType&& c) {
    return internals::unary_operation_modifier<std::decay_t<CallableType>>(
        std::forward<CallableType>(c));
  }  // end of transform

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_UNARYOPERATION_IXX */
