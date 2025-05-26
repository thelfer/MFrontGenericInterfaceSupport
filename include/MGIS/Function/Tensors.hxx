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
  struct FunctionResultTypeTraits<TensorType&> {
    static constexpr auto is_specialized = true;
  };

  template <TensorConcept TensorType>
  struct FunctionResultTypeTraits<tfel::math::View<TensorType>> {
    static constexpr auto is_specialized = true;
  };

  template <TensorConcept TensorType>
  struct tensor_modifier {
    template <FunctionConcept FunctionType>
    auto operator()(FunctionType&) const
        requires(number_of_components<FunctionType> == dynamic_extent
                     ? true
                     : compile_time_size<TensorType> ==
                           number_of_components<FunctionType>);

    template <EvaluatorConcept EvaluatorType>
    auto operator()(const EvaluatorType&) const requires(
        (areTensorModifierRequirementsSatisfied<TensorType,
                                                std::decay_t<EvaluatorType>>));
  };

  struct RotateModifier {
    //
    using CallableType = decltype(
        []<typename TensorType>(const TensorType& t,
                                const tfel::math::tmatrix<3, 3, real>& R)  //
        requires((tfel::math::TensorConcept<TensorType>) ||
                 (tfel::math::StensorConcept<TensorType>)) {
          return tfel::math::change_basis(t, R);
        });

    inline auto operator()(const tfel::math::tmatrix<3, 3, real>& R) const {
      auto c = [R]<typename TensorType>(const TensorType& t)  //
          requires((tfel::math::TensorConcept<TensorType>) ||
                   (tfel::math::StensorConcept<TensorType>)) {
        return tfel::math::change_basis(t, R);
      };
      return internals::unary_operation_modifier<decltype(c)>(c);
    }

    template <typename SecondEvaluatorType>
    auto operator()(SecondEvaluatorType&& e2) const
        requires(EvaluatorConcept<std::decay_t<SecondEvaluatorType>>) {
      return BinaryOperatorCurrying2<CallableType,
                                     std::decay_t<SecondEvaluatorType>>(
          std::forward<SecondEvaluatorType>(e2));
    }  // end of operator()

    template <typename FirstEvaluatorType, typename SecondEvaluatorType>
    auto operator()(FirstEvaluatorType&& e1, SecondEvaluatorType&& e2) const
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
  };

  struct RotateBackwardsModifier {
    //
    using CallableType = decltype(
        []<typename TensorType>(const TensorType& t,
                                const tfel::math::tmatrix<3, 3, real>& R)  //
        requires((tfel::math::TensorConcept<TensorType>) ||
                 (tfel::math::StensorConcept<TensorType>)) {
          return tfel::math::change_basis(t, tfel::math::transpose(R));
        });

    inline auto operator()(const tfel::math::tmatrix<3, 3, real>& R) const {
      auto c = [Rb = tfel::math::transpose(R)]<typename TensorType>(
          const TensorType& t)  //
          requires((tfel::math::TensorConcept<TensorType>) ||
                   (tfel::math::StensorConcept<TensorType>)) {
        return tfel::math::change_basis(t, Rb);
      };
      return internals::unary_operation_modifier<decltype(c)>(c);
    }

    template <typename SecondEvaluatorType>
    auto operator()(SecondEvaluatorType&& e2) const
        requires(EvaluatorConcept<std::decay_t<SecondEvaluatorType>>) {
      return BinaryOperatorCurrying2<CallableType,
                                     std::decay_t<SecondEvaluatorType>>(
          std::forward<SecondEvaluatorType>(e2));
    }  // end of operator()

    template <typename FirstEvaluatorType, typename SecondEvaluatorType>
    auto operator()(FirstEvaluatorType&& e1, SecondEvaluatorType&& e2) const
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
  };   // end of RotateBackwardModifier

}  // end of namespace mgis::function::internals

namespace mgis::function {

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  auto operator|(FunctionType&,
                 const internals::tensor_modifier<TensorType>&)  //
      requires(number_of_components<std::decay_t<FunctionType>> ==
                        dynamic_extent
                    ? true
                    : compile_time_size<TensorType> ==
                          number_of_components<std::decay_t<FunctionType>>);

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
        auto s2 = tfel::math::EvaluationResult<decltype(s)>(s);
        return s2.template computeEigenValues<esolver>();
      });

  inline constexpr auto rotate = internals::RotateModifier{};

  inline constexpr auto rotate_backwards = internals::RotateBackwardsModifier{};

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
