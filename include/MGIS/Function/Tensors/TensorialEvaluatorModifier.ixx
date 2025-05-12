/*!
 * \file   MGIS/Function/TensorialEvaluatorModifier.ixx
 * \brief    
 * \author Thomas Helfer
 * \date   10/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_TENSORIALEVALUATORMODIFIER_IXX
#define LIB_MGIS_FUNCTION_TENSORIALEVALUATORMODIFIER_IXX

#include <utility>

namespace mgis::function::internals {

  template <TensorialObjectConcept TensorType, EvaluatorConcept EvaluatorType>
  constexpr size_type
  TensorialEvaluatorModifier<TensorType,
                             EvaluatorType>::getNumberOfComponents() noexcept {
    return compile_time_size<TensorType>;
  }

  template <TensorialObjectConcept TensorType, EvaluatorConcept EvaluatorType>
  auto TensorialEvaluatorModifier<TensorType, EvaluatorType>::apply(
      const evaluator_result<EvaluatorType>& values) const {
    return tfel::math::map<const TensorType>(values.data());
  };  // end of apply

  template <TensorialObjectConcept TensorType>
  template <typename EvaluatorType>
  auto tensor_modifier<TensorType>::operator()(EvaluatorType&& e) const
      requires((EvaluatorConcept<std::decay_t<EvaluatorType>>)&&(
          areTensorialEvaluatorModifierRequirementsSatisfied<TensorType,
                                                             EvaluatorType>)) {
    using Modifier =
        TensorialEvaluatorModifier<TensorType, std::decay_t<EvaluatorType>>;
    return Modifier(std::forward<EvaluatorType>(e));
  }  // end of operator()

}  // end of namespace mgis::function::internals

#endif /* LIB_MGIS_FUNCTION_TENSORIALEVALUATORMODIFIER_IXX */
