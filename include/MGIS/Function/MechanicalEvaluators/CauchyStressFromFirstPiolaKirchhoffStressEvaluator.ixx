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
  void CauchyStressFromFirstPiolaKirchhoffStressEvaluator<
      N,
      PK1EvaluatorType,
      DeformationGradientEvaluatorType>::check() const {
    using namespace tfel::math;
    BinaryOperationEvaluatorBase<
        CauchyStressFromFirstPiolaKirchhoffStressEvaluator<
            N, PK1EvaluatorType, DeformationGradientEvaluatorType>,
        PK1EvaluatorType, DeformationGradientEvaluatorType>::check();
    const auto nPK1 = this->first_evaluator.getNumberOfComponents();
    raise_if(nPK1 != TensorDimeToSize<N>::value,
             "CauchyStressFromFirstPiolaKirchhoffStressEvaluator::check: "
             "incompatible number of components of the first Piola-Kirchhoff "
             "stress");
    const auto nF = this->second_evaluator.getNumberOfComponents();
    raise_if(nF != TensorDimeToSize<N>::value,
             "CauchyStressFromFirstPiolaKirchhoffStressEvaluator::check: "
             "incompatible number of components of the deformation gradient");
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

}  // end namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_CAUCHYSTRESSFROMFIRSTPIOLAKIRCHHOFFSTRESSEVALUATOR_IXX \
        */
