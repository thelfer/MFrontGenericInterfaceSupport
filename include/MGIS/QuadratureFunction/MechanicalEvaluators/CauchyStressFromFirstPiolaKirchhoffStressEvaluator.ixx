/*!
 * \file   CauchyStressFromFirstPiolaKirchhoffStressEvaluator.ixx
 * \brief    
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#ifndef LIB_MGIS_QUADRATUREFUNCTION_CAUCHYSTRESSFROMFIRSTPIOLAKIRCHHOFFSTRESSEVALUATOR_IXX
#define LIB_MGIS_QUADRATUREFUNCTION_CAUCHYSTRESSFROMFIRSTPIOLAKIRCHHOFFSTRESSEVALUATOR_IXX

#include "TFEL/Math/tensor.hxx"
#include "TFEL/Math/stensor.hxx"
#include "TFEL/Math/Array/View.hxx"
#include "MGIS/Raise.hxx"

namespace mgis::quadrature_function{

  template <unsigned short N,
            EvaluatorConcept DeformationGradientEvaluatorType,
            EvaluatorConcept PK1EvaluatorType>
  CauchyStressFromFirstPiolaKirchhoffStressEvaluator<
      N,
      DeformationGradientEvaluatorType,
      PK1EvaluatorType>::
      CauchyStressFromFirstPiolaKirchhoffStressEvaluator(
          const DeformationGradientEvaluatorType& e1,
          const PK1EvaluatorType& e2)
      : deformation_gradient_evaluator(e1),
        pk1_evaluator(e2) {
  }  // end of CauchyStressFromFirstPiolaKirchhoffStressEvaluator

  template <unsigned short N,
            EvaluatorConcept DeformationGradientEvaluatorType,
            EvaluatorConcept PK1EvaluatorType>
  void CauchyStressFromFirstPiolaKirchhoffStressEvaluator<
      N,
      DeformationGradientEvaluatorType,
      PK1EvaluatorType>::check() const {
    using namespace tfel::math;
    this->deformation_gradient_evaluator.check();
    this->pk1_evaluator.check();
    const auto nF = this->deformation_gradient_evaluator.getNumberOfComponents();
    const auto nPK1 = this->pk1_evaluator.getNumberOfComponents();
    checkMatchingAbstractQuadratureSpaces(this->deformation_gradient_evaluator,
                                          this->pk1_evaluator);
    raise_if(nF != TensorDimeToSize<N>::value,
             "CauchyStressFromFirstPiolaKirchhoffStressEvaluator::check: "
             "incompatible number of components of the deformation gradient");
    raise_if(nPK1 != TensorDimeToSize<N>::value,
             "CauchyStressFromFirstPiolaKirchhoffStressEvaluator::check: "
             "incompatible number of components of the first Piola-Kirchhoff "
             "stress");
  }  // end of check

  template <unsigned short N,
            EvaluatorConcept DeformationGradientEvaluatorType,
            EvaluatorConcept PK1EvaluatorType>
  void CauchyStressFromFirstPiolaKirchhoffStressEvaluator<
      N,
      DeformationGradientEvaluatorType,
      PK1EvaluatorType>::allocateWorkspace() {
    this->deformation_gradient_evaluator.allocatWorkspace();
    this->pk1_evaluator.allocatWorkspace();
  }  // end of allocatWorkspace

  template <unsigned short N,
            EvaluatorConcept DeformationGradientEvaluatorType,
            EvaluatorConcept PK1EvaluatorType>
  const AbstractQuadratureSpace&
  CauchyStressFromFirstPiolaKirchhoffStressEvaluator<
      N,
      DeformationGradientEvaluatorType,
      PK1EvaluatorType>::getQuadratureSpace() const {
    return this->deformation_gradient_evaluator.getQuadratureSpace();
  }  // end of getQuadratureSpace

  template <unsigned short N,
            EvaluatorConcept DeformationGradientEvaluatorType,
            EvaluatorConcept PK1EvaluatorType>
  constexpr size_type CauchyStressFromFirstPiolaKirchhoffStressEvaluator<
      N,
      DeformationGradientEvaluatorType,
      PK1EvaluatorType>::getNumberOfComponents() const noexcept {
    return tfel::math::StensorDimeToSize<N>::value;
  }  // end of getNumberOfComponents

  template <unsigned short N,
            EvaluatorConcept DeformationGradientEvaluatorType,
            EvaluatorConcept PK1EvaluatorType>
  real CauchyStressFromFirstPiolaKirchhoffStressEvaluator<
      N,
      DeformationGradientEvaluatorType,
      PK1EvaluatorType>::operator()(const size_type i) const {
    using namespace tfel::math;
    const auto F =
        map<tensor<N, real>>(this->deformation_gradient_evaluator(i).data());
    const auto pk1 = map<tensor<N, real>>(this->pk1_evaluator(i).data());
    return convertFirstPiolaKirchhoffStressToCauchyStress(pk1, F);
  } // end of operator()

}  // end namespace mgis::quadrature_function

#endif /* LIB_MGIS_QUADRATUREFUNCTION_CAUCHYSTRESSFROMFIRSTPIOLAKIRCHHOFFSTRESSEVALUATOR_IXX */
