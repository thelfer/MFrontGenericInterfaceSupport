/*!
 * \file   MGIS/Function/Tensors.hxx
 * \brief
 * \author Thomas Helfer
 * \date   11/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_TENSORS_HXX
#define LIB_MGIS_FUNCTION_TENSORS_HXX

#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/Evaluator.hxx"

#ifdef MGIS_HAVE_TFEL

#include "MGIS/Function/Tensors/TensorConcept.hxx"
#include "MGIS/Function/Tensors/TensorView.hxx"
#include "MGIS/Function/Tensors/TensorModifier.hxx"

namespace mgis::function::internals {

  template <TensorConcept TensorType>
  struct tensor_modifier {
    template <FunctionalSpaceConcept Space,
              DataLayoutDescription layout,
              bool is_mutable>
    auto operator()(const FunctionView<Space, layout, is_mutable>&) const
        requires(layout.data_size == dynamic_extent
                     ? true
                     : compile_time_size<TensorType> == layout.data_size);

    template <FunctionalSpaceConcept Space, size_type N>
    auto operator()(const Function<Space, N>&) const
        requires(N == dynamic_extent ? true
                                     : compile_time_size<TensorType> == N);

    template <FunctionalSpaceConcept Space, size_type N>
    auto operator()(Function<Space, N>&) const
        requires(N == dynamic_extent ? true
                                     : compile_time_size<TensorType> == N);

    template <typename EvaluatorType>
    auto operator()(EvaluatorType&&) const
        requires((EvaluatorConcept<std::decay_t<EvaluatorType>>)&&(
            areTensorModifierRequirementsSatisfied<TensorType, EvaluatorType>));
  };

}  // end of namespace mgis::function::internals

namespace mgis::function {

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable,
            TensorConcept TensorType>
  auto operator|(const FunctionView<Space, layout, is_mutable>&,
                 const internals::tensor_modifier<TensorType>&)  //
      requires(layout.data_size == dynamic_extent
                   ? true
                   : compile_time_size<TensorType> == layout.data_size);

  template <FunctionalSpaceConcept Space, size_type N, TensorConcept TensorType>
  auto operator|(const Function<Space, N>&,
                 const internals::tensor_modifier<TensorType>&)  //
      requires(N == dynamic_extent ? true : compile_time_size<TensorType> == N);

  template <FunctionalSpaceConcept Space, size_type N, TensorConcept TensorType>
  auto operator|(Function<Space, N>& f,
                 const internals::tensor_modifier<TensorType>&)  //
      requires(N == dynamic_extent ? true : compile_time_size<TensorType> == N);

  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      inline constexpr auto as_stensor =
          internals::tensor_modifier<tfel::math::stensor<N, real>>{};

  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      inline constexpr auto as_tensor =
          internals::tensor_modifier<tfel::math::tensor<N, real>>{};

  template <unsigned short N>
  inline constexpr auto as_fsarray =
      internals::tensor_modifier<tfel::math::fsarray<N, real>>{};

  template <unsigned short N>
  inline constexpr auto as_tvector =
      internals::tensor_modifier<tfel::math::tvector<N, real>>{};

  template <unsigned short N, unsigned short M>
  inline constexpr auto as_tmatrix =
      internals::tensor_modifier<tfel::math::tmatrix<N, M, real>>{};

  inline constexpr auto trace = internals::unary_operation_modifier2(
      []<typename TensorType>(const TensorType& t) requires(
          (tfel::math::TensorConcept<TensorType>) ||
          (tfel::math::StensorConcept<TensorType>)) {
        return tfel::math::trace(t);
      });

  inline constexpr auto det = internals::unary_operation_modifier2(
      []<typename TensorType>(const TensorType& t) requires(
          (tfel::math::TensorConcept<TensorType>) ||
          (tfel::math::StensorConcept<TensorType>)) {
        return tfel::math::det(t);
      });

  template <tfel::math::stensor_common::EigenSolver esolver =
                tfel::math::stensor_common::TFELEIGENSOLVER>
  inline constexpr auto eigen_values = internals::unary_operation_modifier2(
      [](const tfel::math::StensorConcept auto& s) {
        auto s2 = eval(s);
        return s2.template computeEigenValues<esolver>();
      });

  inline constexpr auto rotate = internals::binary_operation_modifier2(
      []<typename TensorType>(const TensorType& pk1,
                              const tfel::math::tmatrix<3, 3>& R)  //
      requires((tfel::math::TensorConcept<TensorType>) ||
               (tfel::math::StensorConcept<TensorType>)) {
        return tfel::math::change_basis(pk1, R);
      });

  inline constexpr auto rotate_backwards =
      internals::binary_operation_modifier2(
          []<typename TensorType>(const TensorType& pk1,
                                  const tfel::math::tmatrix<3, 3>& R)  //
          requires((tfel::math::TensorConcept<TensorType>) ||
                   (tfel::math::StensorConcept<TensorType>)) {
            return tfel::math::change_basis(pk1, transpose(R));
          });

  template <tfel::math::stensor_common::EigenSolver esolver =
                tfel::math::stensor_common::TFELEIGENSOLVER>
  inline constexpr auto logarithmic_strain =
      internals::unary_operation_modifier2(
          [](const tfel::math::TensorConcept auto& F) {
            using namespace tfel::math;
            constexpr auto N = getSpaceDimension<decltype(F)>();
            const auto e = computeGreenLagrangeTensor(F);
            const auto [vp, m] = e.template computeEigenVectors<esolver>();
            const auto log_vp =
                map([](const real x) { return std::log1p(2 * x) / 2; }, vp);
            return stensor<N, real>::computeIsotropicFunction(e, m);
          });

}  // end of namespace mgis::function

#endif MGIS_HAVE_TFEL

#include "MGIS/Function/Tensors.ixx"

#endif /* LIB_MGIS_FUNCTION_TENSORS_HXX */
