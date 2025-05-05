/*!
 * \file   HydrostaticStressEvaluator.ixx
 * \brief
 * \author Thomas Helfer
 * \date   01/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_HYDROSTATICSTRESSSTRESSEVALUATOR_IXX
#define LIB_MGIS_FUNCTION_HYDROSTATICSTRESSSTRESSEVALUATOR_IXX

#include "TFEL/Math/stensor.hxx"
#include "TFEL/Math/Array/View.hxx"

namespace mgis::function {

  template <unsigned short N, EvaluatorConcept StressEvaluatorType>
  constexpr size_type
  HydrostaticStressEvaluator<N, StressEvaluatorType>::getNumberOfComponents()
      const noexcept {
    return 1;
  }  // end of getNumberOfComponents

  template <unsigned short N, EvaluatorConcept StressEvaluatorType>
  real HydrostaticStressEvaluator<N, StressEvaluatorType>::operator()(
      const size_type i) const {
    using namespace tfel::math;
    const auto sig = map<stensor<N, real>>(this->stress_evaluator(i).data());
    return trace(sig) / 3;
  }

  template <unsigned short N, EvaluatorConcept StressEvaluatorType>
  auto hydrostatic_stress(const StressEvaluatorType& e) requires((N == 1) ||
                                                                 (N == 2) ||
                                                                 (N == 3)) {
    return HydrostaticStressEvaluator<N, StressEvaluatorType>(e);
  }  // end of hydrostatic_stress

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_HYDROSTATICSTRESSSTRESSEVALUATOR_IXX */
