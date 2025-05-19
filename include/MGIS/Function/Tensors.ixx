/*!
 * \file   MGIS/Function/Tensors.ixx
 * \brief
 * \author Thomas Helfer
 * \date   13/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_TENSORS_IXX
#define LIB_MGIS_FUNCTION_TENSORS_IXX

#include <iostream>

namespace mgis::function::internals {

  template <TensorConcept TensorType>
  template <FunctionConcept FunctionType>
  auto tensor_modifier<TensorType>::operator()(FunctionType& f) const
      requires(number_of_components<FunctionType> == dynamic_extent
                   ? true
                   : compile_time_size<TensorType> ==
                         number_of_components<FunctionType>) {
    return TensorView<FunctionType, TensorType>(f);
  }

  template <TensorConcept TensorType>
  template <typename EvaluatorType>
  auto tensor_modifier<TensorType>::operator()(EvaluatorType&& e) const
      requires((EvaluatorConcept<std::decay_t<EvaluatorType>>)&&(
          areTensorModifierRequirementsSatisfied<
              TensorType,
              std::decay_t<EvaluatorType>>)) {
    using Modifier = TensorModifier<TensorType, std::decay_t<EvaluatorType>>;
    return Modifier(std::forward<EvaluatorType>(e));
  }  // end of operator()

}  // end of namespace mgis::function::internals

namespace mgis::function {

  template <typename FunctionType, TensorConcept TensorType>
  auto operator|(FunctionType&& f,
                 const internals::tensor_modifier<TensorType>& m)  //
      requires((FunctionConcept<std::decay_t<FunctionType>>)&&     //
               (!std::is_rvalue_reference_v<FunctionType&&>)&&     //
               (number_of_components<std::decay_t<FunctionType>> ==
                        dynamic_extent
                    ? true
                    : compile_time_size<TensorType> ==
                          number_of_components<std::decay_t<FunctionType>>)) {
    return m(f);
  }

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TENSORS_IXX */
