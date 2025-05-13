/*!
 * \file   MGIS/Function/TensorModifier.ixx
 * \brief    
 * \author Thomas Helfer
 * \date   13/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_TENSORMODIFIER_IXX
#define LIB_MGIS_FUNCTION_TENSORMODIFIER_IXX

namespace mgis::function {

  template <TensorConcept TensorType, EvaluatorConcept EvaluatorType>
  requires(areTensorModifierRequirementsSatisfied<TensorType,
                                                  EvaluatorType>)  //
      constexpr size_type
      TensorModifier<TensorType,
                     EvaluatorType>::getNumberOfComponents() noexcept {
    return compile_time_size<TensorType>;
  }

  template <TensorConcept TensorType, EvaluatorConcept EvaluatorType>
  requires(areTensorModifierRequirementsSatisfied<TensorType,
                                                  EvaluatorType>)  //
      auto TensorModifier<TensorType, EvaluatorType>::apply(
          const evaluator_result<EvaluatorType>& values) const {
    return tfel::math::map<const TensorType>(values.data());
  }  // end of apply

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TENSORMODIFIER_IXX */
