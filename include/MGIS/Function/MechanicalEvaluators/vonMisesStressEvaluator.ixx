/*!
 * \file   vonMisesStressEvaluator.ixx
 * \brief
 * \author Thomas Helfer
 * \date   01/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_VONMISESSTRESSEVALUATOR_IXX
#define LIB_MGIS_FUNCTION_VONMISESSTRESSEVALUATOR_IXX

#include <utility>
#include "TFEL/Math/stensor.hxx"
#include "TFEL/Math/Array/View.hxx"

namespace mgis::function {

  template <unsigned short N, EvaluatorConcept StressEvaluatorType>
  constexpr size_type
  vonMisesStressEvaluator<N, StressEvaluatorType>::getNumberOfComponents()
      const noexcept {
    return 1;
  }  // end of getNumberOfComponents

  template <unsigned short N, EvaluatorConcept StressEvaluatorType>
  real vonMisesStressEvaluator<N, StressEvaluatorType>::operator()(
      const size_type i) const {
    using namespace tfel::math;
    const auto sig = map<stensor<N, real>>(this->stress_evaluator(i).data());
    return sigmaeq(sig);
  }

  template <unsigned short N>
  template <EvaluatorConcept StressEvaluatorType>
  constexpr auto vmis_fn<N>::operator()(StressEvaluatorType&& e) const {
    return vonMisesStressEvaluator<N, StressEvaluatorType>(
        std::forward<StressEvaluatorType>(e));
  }  // end of vmis

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_VONMISESSTRESSEVALUATOR_IXX */
