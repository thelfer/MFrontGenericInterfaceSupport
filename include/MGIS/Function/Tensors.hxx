/*!
 * \file   MGIS/Function/Tensors.hxx
 * \brief
 * \author Thomas Helfer
 * \date   11/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_TENSORS_HXX
#define LIB_MGIS_FUNCTION_TENSORS_HXX

#include "MGIS/Function/FunctionConcept.hxx"
#include "MGIS/Function/EvaluatorConcept.hxx"

#ifdef MGIS_HAVE_TFEL

#include "MGIS/Function/Tensors/TensorConcept.hxx"
#include "MGIS/Function/Tensors/TensorView.hxx"
#include "MGIS/Function/Tensors/TensorModifier.hxx"
#include "MGIS/Function/Tensors/CoalescedMemoryAccessTensorView.hxx"
#include "MGIS/Function/Tensors/CoalescedMemoryAccessCompositeTensorsView.hxx"
#include "MGIS/Function/Tensors/StridedCoalescedMemoryAccessTensorView.hxx"
#include "MGIS/Function/Tensors/StridedCoalescedMemoryAccessCompositeTensorsView.hxx"

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
    //! \brief this alias allows to match the Evaluator Modifier concept
    using Tag = ::mgis::function::EvaluatorModifierTag;

    template <FunctionConcept FunctionType>
    constexpr auto operator()(FunctionType&) const
        requires(number_of_components<FunctionType> == dynamic_extent
                     ? true
                     : compile_time_size<TensorType> ==
                           number_of_components<FunctionType>);
    /*!
     * \brief create a new modifier
     * \param[in] e: evaluator type
     */
    template <EvaluatorConcept EvaluatorType>
    constexpr auto operator()(const EvaluatorType&) const requires(
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

    constexpr auto operator()(const tfel::math::tmatrix<3, 3, real>& R) const {
      auto c = [R]<typename TensorType>(const TensorType& t)  //
          requires((tfel::math::TensorConcept<TensorType>) ||
                   (tfel::math::StensorConcept<TensorType>)) {
        return tfel::math::change_basis(t, R);
      };
      return internals::unary_operation_modifier<decltype(c)>(c);
    }

    template <typename SecondEvaluatorType>
    constexpr auto operator()(SecondEvaluatorType&& e2) const
        requires(EvaluatorConcept<std::decay_t<SecondEvaluatorType>>) {
      return BinaryOperatorCurrying2<CallableType,
                                     std::decay_t<SecondEvaluatorType>>(
          std::forward<SecondEvaluatorType>(e2));
    }  // end of operator()

    template <typename FirstEvaluatorType, typename SecondEvaluatorType>
    constexpr auto operator()(FirstEvaluatorType&& e1,
                              SecondEvaluatorType&& e2) const
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

    constexpr auto operator()(const tfel::math::tmatrix<3, 3, real>& R) const {
      auto c = [Rb = tfel::math::transpose(R)]<typename TensorType>(
          const TensorType& t)  //
          requires((tfel::math::TensorConcept<TensorType>) ||
                   (tfel::math::StensorConcept<TensorType>)) {
        return tfel::math::change_basis(t, Rb);
      };
      return internals::unary_operation_modifier<decltype(c)>(c);
    }

    template <typename SecondEvaluatorType>
    constexpr auto operator()(SecondEvaluatorType&& e2) const
        requires(EvaluatorConcept<std::decay_t<SecondEvaluatorType>>) {
      return BinaryOperatorCurrying2<CallableType,
                                     std::decay_t<SecondEvaluatorType>>(
          std::forward<SecondEvaluatorType>(e2));
    }  // end of operator()

    template <typename FirstEvaluatorType, typename SecondEvaluatorType>
    constexpr auto operator()(FirstEvaluatorType&& e1,
                              SecondEvaluatorType&& e2) const
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

  template <typename EvaluatorType>
  concept Stensor1DEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&(
      std::is_convertible_v<evaluator_result<EvaluatorType>,
                            tfel::math::stensor<1u, real>>);

  template <typename EvaluatorType>
  concept Stensor2DEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&(
      std::is_convertible_v<evaluator_result<EvaluatorType>,
                            tfel::math::stensor<2u, real>>);

  template <typename EvaluatorType>
  concept Stensor3DEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&(
      std::is_convertible_v<evaluator_result<EvaluatorType>,
                            tfel::math::stensor<3u, real>>);

  template <typename EvaluatorType>
  concept StensorEvaluatorConcept =
      (Stensor1DEvaluatorConcept<EvaluatorType>) ||
      (Stensor2DEvaluatorConcept<EvaluatorType>) ||
      (Stensor3DEvaluatorConcept<EvaluatorType>);

  template <typename EvaluatorType>
  concept Tensor1DEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&(
      std::is_convertible_v<evaluator_result<EvaluatorType>,
                            tfel::math::tensor<1u, real>>);

  template <typename EvaluatorType>
  concept Tensor2DEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&(
      std::is_convertible_v<evaluator_result<EvaluatorType>,
                            tfel::math::tensor<2u, real>>);

  template <typename EvaluatorType>
  concept Tensor3DEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&(
      std::is_convertible_v<evaluator_result<EvaluatorType>,
                            tfel::math::tensor<3u, real>>);

  template <typename EvaluatorType>
  concept TensorEvaluatorConcept = (Tensor1DEvaluatorConcept<EvaluatorType>) ||
                                   (Tensor2DEvaluatorConcept<EvaluatorType>) ||
                                   (Tensor3DEvaluatorConcept<EvaluatorType>);

  template <typename EvaluatorType>
  concept ST2toST21DEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&(
      std::is_convertible_v<evaluator_result<EvaluatorType>,
                            tfel::math::st2tost2<1u, real>>);

  template <typename EvaluatorType>
  concept ST2toST22DEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&(
      std::is_convertible_v<evaluator_result<EvaluatorType>,
                            tfel::math::st2tost2<2u, real>>);

  template <typename EvaluatorType>
  concept ST2toST23DEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&(
      std::is_convertible_v<evaluator_result<EvaluatorType>,
                            tfel::math::st2tost2<3u, real>>);

  template <typename EvaluatorType>
  concept ST2toST2EvaluatorConcept =
      (ST2toST21DEvaluatorConcept<EvaluatorType>) ||
      (ST2toST22DEvaluatorConcept<EvaluatorType>) ||
      (ST2toST23DEvaluatorConcept<EvaluatorType>);

  template <typename EvaluatorType>
  concept ST2toT21DEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&(
      std::is_convertible_v<evaluator_result<EvaluatorType>,
                            tfel::math::st2tot2<1u, real>>);

  template <typename EvaluatorType>
  concept ST2toT22DEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&(
      std::is_convertible_v<evaluator_result<EvaluatorType>,
                            tfel::math::st2tot2<2u, real>>);

  template <typename EvaluatorType>
  concept ST2toT23DEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&(
      std::is_convertible_v<evaluator_result<EvaluatorType>,
                            tfel::math::st2tot2<3u, real>>);

  template <typename EvaluatorType>
  concept ST2toT2EvaluatorConcept =
      (ST2toT21DEvaluatorConcept<EvaluatorType>) ||
      (ST2toT22DEvaluatorConcept<EvaluatorType>) ||
      (ST2toT23DEvaluatorConcept<EvaluatorType>);

  template <typename EvaluatorType>
  concept T2toST21DEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&(
      std::is_convertible_v<evaluator_result<EvaluatorType>,
                            tfel::math::t2tost2<1u, real>>);

  template <typename EvaluatorType>
  concept T2toST22DEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&(
      std::is_convertible_v<evaluator_result<EvaluatorType>,
                            tfel::math::t2tost2<2u, real>>);

  template <typename EvaluatorType>
  concept T2toST23DEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&(
      std::is_convertible_v<evaluator_result<EvaluatorType>,
                            tfel::math::t2tost2<3u, real>>);

  template <typename EvaluatorType>
  concept T2toST2EvaluatorConcept =
      (T2toST21DEvaluatorConcept<EvaluatorType>) ||
      (T2toST22DEvaluatorConcept<EvaluatorType>) ||
      (T2toST23DEvaluatorConcept<EvaluatorType>);

  template <typename EvaluatorType>
  concept T2toT21DEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&(
      std::is_convertible_v<evaluator_result<EvaluatorType>,
                            tfel::math::t2tot2<1u, real>>);

  template <typename EvaluatorType>
  concept T2toT22DEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&(
      std::is_convertible_v<evaluator_result<EvaluatorType>,
                            tfel::math::t2tot2<2u, real>>);

  template <typename EvaluatorType>
  concept T2toT23DEvaluatorConcept = (EvaluatorConcept<EvaluatorType>)&&(
      std::is_convertible_v<evaluator_result<EvaluatorType>,
                            tfel::math::t2tot2<3u, real>>);

  template <typename EvaluatorType>
  concept T2toT2EvaluatorConcept = (T2toT21DEvaluatorConcept<EvaluatorType>) ||
                                   (T2toT22DEvaluatorConcept<EvaluatorType>) ||
                                   (T2toT23DEvaluatorConcept<EvaluatorType>);

  template <typename EvaluatorType>
  concept FourthOrderTensorEvaluatorConcept =
      ST2toST2EvaluatorConcept<EvaluatorType> ||
      T2toST2EvaluatorConcept<EvaluatorType> ||
      ST2toT2EvaluatorConcept<EvaluatorType> ||
      T2toT2EvaluatorConcept<EvaluatorType>;

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr auto operator|(FunctionType&,
                           const internals::tensor_modifier<TensorType>&)  //
      requires(number_of_components<FunctionType> == dynamic_extent
                   ? true
                   : compile_time_size<TensorType> ==
                         number_of_components<FunctionType>);

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

  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      inline constexpr auto as_st2tost2 =
          internals::tensor_modifier<tfel::math::st2tost2<N, real>>{};

  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      inline constexpr auto as_t2tost2 =
          internals::tensor_modifier<tfel::math::t2tost2<N, real>>{};

  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      inline constexpr auto as_st2tot2 =
          internals::tensor_modifier<tfel::math::st2tot2<N, real>>{};

  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      inline constexpr auto as_t2tot2 =
          internals::tensor_modifier<tfel::math::t2tot2<N, real>>{};

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

#endif /* MGIS_HAVE_TFEL */

#include "MGIS/Function/Tensors.ixx"

#endif /* LIB_MGIS_FUNCTION_TENSORS_HXX */
