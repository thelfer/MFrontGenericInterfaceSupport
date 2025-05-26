/*!
 * \file   MGIS/Function/Tensors.ixx
 * \brief
 * \author Thomas Helfer
 * \date   13/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_TENSORS_IXX
#define LIB_MGIS_FUNCTION_TENSORS_IXX

namespace mgis::function::customization_points {

  template <TensorConcept TensorType>
  struct AbsoluteValue<TensorType> {
    static constexpr auto exe(const TensorType& v) noexcept {
      constexpr auto N = compile_time_size<TensorType>;
      return absolue_value(std::span<const real, N>(v.data(), N));
    }
  };

  template <TensorConcept TensorType>
  struct AbsoluteValue<tfel::math::View<const TensorType>> {
    static constexpr auto exe(
        const tfel::math::View<const TensorType>& v) noexcept {
      constexpr auto N = compile_time_size<TensorType>;
      return absolue_value(std::span<const real, N>(v.data(), N));
    }
  };

  template <TensorConcept TensorType>
  struct MaximumComponent<TensorType> {
    static constexpr real exe(const TensorType& v) noexcept {
      constexpr auto N = compile_time_size<TensorType>;
      return maximum_component(std::span<const real, N>(v.data(), N));
    }
  };

  template <TensorConcept TensorType>
  struct MaximumComponent<tfel::math::View<const TensorType>> {
    static constexpr real exe(
        const tfel::math::View<const TensorType>& v) noexcept {
      constexpr auto N = compile_time_size<TensorType>;
      return maximum_component(std::span<const real, N>(v.data(), N));
    }
  };

  template <TensorConcept TensorType>
  struct MinimumComponent<TensorType> {
    static constexpr real exe(const TensorType& v) noexcept {
      constexpr auto N = compile_time_size<TensorType>;
      return minimum_component(std::span<const real, N>(v.data(), N));
    }
  };

  template <TensorConcept TensorType>
  struct MinimumComponent<tfel::math::View<const TensorType>> {
    static constexpr real exe(
        const tfel::math::View<const TensorType>& v) noexcept {
      constexpr auto N = compile_time_size<TensorType>;
      return minimum_component(std::span<const real, N>(v.data(), N));
    }
  };

}  // end of namespace mgis::function::customization_points

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
  template <EvaluatorConcept EvaluatorType>
  auto tensor_modifier<TensorType>::operator()(const EvaluatorType& e) const
      requires((areTensorModifierRequirementsSatisfied<
                   TensorType,
                   std::decay_t<EvaluatorType>>)) {
    using Modifier = TensorModifier<TensorType, std::decay_t<EvaluatorType>>;
    return Modifier(e);
  }  // end of operator()

}  // end of namespace mgis::function::internals

namespace mgis::function {

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  auto operator|(FunctionType& f,
                 const internals::tensor_modifier<TensorType>& m)  //
      requires(number_of_components<std::decay_t<FunctionType>> ==
                        dynamic_extent
                    ? true
                    : compile_time_size<TensorType> ==
                          number_of_components<std::decay_t<FunctionType>>) {
    return m(f);
  }

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TENSORS_IXX */
