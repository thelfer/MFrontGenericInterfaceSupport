/*!
 * \file   MGIS/Function/TransformEvaluatorModifier.ixx
 * \brief
 * \author Thomas Helfer
 * \date   09/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_TRANSFORMEVALUATORMODIFIER_IXX
#define LIB_MGIS_FUNCTION_TRANSFORMEVALUATORMODIFIER_IXX

namespace mgis::function {

  template <typename CallableType, EvaluatorConcept EvaluatorType>
  TransformEvaluatorModifier<CallableType, EvaluatorType>::
      TransformEvaluatorModifier(const CallableType& c, const EvaluatorType& e)
      : EvaluatorModifierBase<
            TransformEvaluatorModifier<CallableType, EvaluatorType>,
            EvaluatorType>(e),
        modifier(c) {}  // end of TransformEvaluatorModifier

  template <typename CallableType, EvaluatorConcept EvaluatorType>
  auto TransformEvaluatorModifier<CallableType, EvaluatorType>::apply(
      const evaluator_result<EvaluatorType>& values) const {
    return this->modifier(values);
  }  // end of apply

  template <typename CallableType, EvaluatorConcept EvaluatorType>
  auto TransformEvaluatorModifier2<CallableType, EvaluatorType>::apply(
      const evaluator_result<EvaluatorType>& values) const {
    auto c = CallableType{};
    return c(values);
  }  // end of apply

  namespace internals {

    template <typename CallableType>
    transform_modifier<CallableType>::transform_modifier(CallableType&& c)
        : modifier(std::forward<CallableType>(c)) {}

    template <typename CallableType>
    template <typename EvaluatorType>
    auto transform_modifier<CallableType>::operator()(EvaluatorType&& e) const
        requires((EvaluatorConcept<std::decay_t<EvaluatorType>>)&&(
            std::invocable<CallableType,
                           evaluator_result<std::decay_t<EvaluatorType>>>)) {
      return TransformEvaluatorModifier<CallableType,
                                        std::decay_t<EvaluatorType>>(
          this->modifier, std::forward<EvaluatorType>(e));
    }  // end of operator()

    template <typename CallableType>
    template <typename EvaluatorType>
    auto transform_modifier2_impl<CallableType>::operator()(EvaluatorType&& e) const
        requires((EvaluatorConcept<std::decay_t<EvaluatorType>>)&&(
            std::invocable<CallableType,
                           evaluator_result<std::decay_t<EvaluatorType>>>)) {
      return TransformEvaluatorModifier2<CallableType,
                                         std::decay_t<EvaluatorType>>(
          std::forward<EvaluatorType>(e));
    }  // end of operator()

    template <typename CallableType>
    constexpr auto transform_modifier2(CallableType){
      return transform_modifier2_impl<CallableType>{};
    } // end of transform_modifier2

  }  // namespace internals

  template <typename CallableType>
  auto transform(CallableType&& c) {
    return internals::transform_modifier<std::decay_t<CallableType>>(
        std::forward<CallableType>(c));
  }  // end of transform

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TRANSFORMEVALUATORMODIFIER_IXX */
