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
  real HydrostaticStressEvaluator<N, StressEvaluatorType>::apply(
      const auto& values) const {
    using namespace tfel::math;
    const auto sig = map<stensor<N, real>>(values.data());
    return trace(sig) / 3;
  }

  template <unsigned short N, typename StressEvaluatorType>
  auto hydrostatic_stress(StressEvaluatorType&& e) requires(
      EvaluatorConcept<std::decay_t<StressEvaluatorType>> &&
      ((N == 1) || (N == 2) || (N == 3))) {
    return HydrostaticStressEvaluator<N, std::decay_t<StressEvaluatorType>>(e);
  }  // end of hydrostatic_stress

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_HYDROSTATICSTRESSSTRESSEVALUATOR_IXX */
