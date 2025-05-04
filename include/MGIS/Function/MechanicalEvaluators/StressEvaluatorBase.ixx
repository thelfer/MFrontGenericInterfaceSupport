/*!
 * \file   StressEvaluatorBase.ixx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_STRESSEVALUATORBASE_IXX
#define LIB_MGIS_FUNCTION_STRESSEVALUATORBASE_IXX

#include "TFEL/Math/tensor.hxx"
#include "TFEL/Math/stensor.hxx"

namespace mgis::function {

  template <unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  StressEvaluatorBase<N, StressEvaluatorType, isSymmetric>::StressEvaluatorBase(
      const StressEvaluatorType& e)
      : stress_evaluator(e) {}  // end of StressEvaluatorBase

  template <unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  StressEvaluatorBase<N, StressEvaluatorType, isSymmetric>::StressEvaluatorBase(
      const StressEvaluatorBase&) = default;

  template <unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  StressEvaluatorBase<N, StressEvaluatorType, isSymmetric>::StressEvaluatorBase(
      StressEvaluatorBase&&) = default;

  template <unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  void StressEvaluatorBase<N, StressEvaluatorType, isSymmetric>::check() const {
    using namespace tfel::math;
    this->stress_evaluator.check();
    const auto nc = this->stress_evaluator.getNumberOfComponents();
    if constexpr (isSymmetric) {
      raise_if(nc != StensorDimeToSize<N>::value,
               "StressEvaluatorBase::check: "
               "incompatible number of components of the stress");
    } else {
      raise_if(nc != TensorDimeToSize<N>::value,
               "StressEvaluatorBase::check: "
               "incompatible number of components of the stress");
    }
  }  // end of check

  template <unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  void StressEvaluatorBase<N, StressEvaluatorType, isSymmetric>::
      allocateWorkspace() {
    this->stress_evaluator.allocatWorkspace();
  }  // end of allocatWorkspace

  template <unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  const auto&
  StressEvaluatorBase<N, StressEvaluatorType, isSymmetric>::getSpace() const {
    return this->stress_evaluator.getSpace();
  }  // end of getSpace

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_STRESSEVALUATORBASE_IXX */
