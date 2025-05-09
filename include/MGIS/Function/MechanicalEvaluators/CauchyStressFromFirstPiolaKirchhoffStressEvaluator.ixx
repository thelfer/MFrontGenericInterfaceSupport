/*!
 * \file   CauchyStressFromFirstPiolaKirchhoffStressEvaluator.ixx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_CAUCHYSTRESSFROMFIRSTPIOLAKIRCHHOFFSTRESSEVALUATOR_IXX
#define LIB_MGIS_FUNCTION_CAUCHYSTRESSFROMFIRSTPIOLAKIRCHHOFFSTRESSEVALUATOR_IXX

#include "TFEL/Math/tensor.hxx"
#include "TFEL/Math/stensor.hxx"
#include "TFEL/Math/Array/View.hxx"
#include "MGIS/Raise.hxx"

namespace mgis::function {

  template <unsigned short N,
            EvaluatorConcept PK1EvaluatorType,
            EvaluatorConcept DeformationGradientEvaluatorType>
  bool CauchyStressFromFirstPiolaKirchhoffStressEvaluator<
      N,
      PK1EvaluatorType,
      DeformationGradientEvaluatorType>::check(Context& ctx) const {
    using namespace tfel::math;
    if (!BinaryOperationEvaluatorBase<
            CauchyStressFromFirstPiolaKirchhoffStressEvaluator<
                N, PK1EvaluatorType, DeformationGradientEvaluatorType>,
            PK1EvaluatorType, DeformationGradientEvaluatorType>::check(ctx)) {
      return false;
    }
    const auto nPK1 = this->first_evaluator.getNumberOfComponents();
    if (nPK1 != TensorDimeToSize<N>::value) {
      return ctx.registerErrorMessage(
          "incompatible number of components of the first Piola-Kirchhoff "
          "stress");
    }
    const auto nF = this->second_evaluator.getNumberOfComponents();
    if (nF != TensorDimeToSize<N>::value) {
      return ctx.registerErrorMessage(
          "incompatible number of components of the deformation gradient");
    }
    return true;
  }  // end of check

  template <unsigned short N,
            EvaluatorConcept PK1EvaluatorType,
            EvaluatorConcept DeformationGradientEvaluatorType>
  constexpr size_type CauchyStressFromFirstPiolaKirchhoffStressEvaluator<
      N,
      PK1EvaluatorType,
      DeformationGradientEvaluatorType>::getNumberOfComponents()
      const noexcept {
    return tfel::math::StensorDimeToSize<N>::value;
  }  // end of getNumberOfComponents

  template <unsigned short N,
            EvaluatorConcept PK1EvaluatorType,
            EvaluatorConcept DeformationGradientEvaluatorType>
  auto CauchyStressFromFirstPiolaKirchhoffStressEvaluator<
      N,
      PK1EvaluatorType,
      DeformationGradientEvaluatorType>::apply(const auto& pk1_values,
                                               const auto& F_values) const {
    using namespace tfel::math;
    const auto F = map<tensor<N, real>>(F_values.data());
    const auto pk1 = map<tensor<N, real>>(pk1_values.data());
    return convertFirstPiolaKirchhoffStressToCauchyStress(pk1, F);
  }  // end of operator()

  namespace internals {

    template <unsigned short N,
              EvaluatorConcept DeformationGradientEvaluatorType>
    from_pk1_to_cauchy_modifier<N, DeformationGradientEvaluatorType>::
        from_pk1_to_cauchy_modifier(DeformationGradientEvaluatorType&& e)
        : F(std::forward<DeformationGradientEvaluatorType>(e)) {}

    template <unsigned short N,
              EvaluatorConcept DeformationGradientEvaluatorType>
    template <typename PK1EvaluatorType>
    auto from_pk1_to_cauchy_modifier<N, DeformationGradientEvaluatorType>::
    operator()(PK1EvaluatorType&& s) const
        requires(EvaluatorConcept<std::decay_t<PK1EvaluatorType>>) {
      return CauchyStressFromFirstPiolaKirchhoffStressEvaluator<
          N, std::decay_t<PK1EvaluatorType>, DeformationGradientEvaluatorType>(
          s, this->F);
    }  // end of operator()

  }  // namespace internals

  template <unsigned short N, typename DeformationGradientEvaluatorType>
  auto from_pk1_to_cauchy(DeformationGradientEvaluatorType&& F) requires(
      EvaluatorConcept<std::decay_t<DeformationGradientEvaluatorType>>) {
    return internals::from_pk1_to_cauchy_modifier<
        N, std::decay_t<DeformationGradientEvaluatorType>>(
        std::forward<DeformationGradientEvaluatorType>(F));
  }

}  // end namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_CAUCHYSTRESSFROMFIRSTPIOLAKIRCHHOFFSTRESSEVALUATOR_IXX \
        */
