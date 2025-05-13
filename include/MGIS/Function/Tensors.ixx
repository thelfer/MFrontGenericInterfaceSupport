/*!
 * \file   MGIS/Function/Tensors.ixx
 * \brief
 * \author Thomas Helfer
 * \date   13/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_TENSORS_IXX
#define LIB_MGIS_FUNCTION_TENSORS_IXX

namespace mgis::function::internals {

  template <TensorConcept TensorType>
  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  auto tensor_modifier<TensorType>::operator()(
      const FunctionView<Space, layout, is_mutable>& f) const
      requires(layout.data_size == dynamic_extent
                   ? true
                   : compile_time_size<TensorType> == layout.data_size) {
    return TensorView<Space, TensorType, is_mutable>(f);
  }

  template <TensorConcept TensorType>
  template <FunctionalSpaceConcept Space, size_type N>
  auto tensor_modifier<TensorType>::operator()(const Function<Space, N>& f)
      const requires(N == dynamic_extent ? true
                                         : compile_time_size<TensorType> == N) {
    return TensorView<Space, TensorType, false>(f.view());
  }

  template <TensorConcept TensorType>
  template <FunctionalSpaceConcept Space, size_type N>
  auto tensor_modifier<TensorType>::operator()(Function<Space, N>& f) const
      requires(N == dynamic_extent ? true
                                   : compile_time_size<TensorType> == N) {
    return TensorView<Space, TensorType, true>(f.view());
  }

  template <TensorConcept TensorType>
  template <typename EvaluatorType>
  auto tensor_modifier<TensorType>::operator()(EvaluatorType&& e) const
      requires((EvaluatorConcept<std::decay_t<EvaluatorType>>)&&(
          areTensorModifierRequirementsSatisfied<TensorType, EvaluatorType>)) {
    using Modifier = TensorModifier<TensorType, std::decay_t<EvaluatorType>>;
    return Modifier(std::forward<EvaluatorType>(e));
  }  // end of operator()

}  // end of namespace mgis::function::internals

namespace mgis::function {

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable,
            TensorConcept TensorType>
  auto operator|(const FunctionView<Space, layout, is_mutable>& f,
                 const internals::tensor_modifier<TensorType>& m)  //
      requires(layout.data_size == dynamic_extent
                   ? true
                   : compile_time_size<TensorType> == layout.data_size) {
    return m(f);
  }

  template <FunctionalSpaceConcept Space, size_type N, TensorConcept TensorType>
  auto operator|(const Function<Space, N>& f,
                 const internals::tensor_modifier<TensorType>& m)  //
      requires(N == dynamic_extent ? true
                                   : compile_time_size<TensorType> == N) {
    return m(f);
  }

  template <FunctionalSpaceConcept Space, size_type N, TensorConcept TensorType>
  auto operator|(Function<Space, N>& f,
                 const internals::tensor_modifier<TensorType>& m)  //
      requires(N == dynamic_extent ? true
                                   : compile_time_size<TensorType> == N) {
    return m(f);
  }

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TENSORS_IXX */
