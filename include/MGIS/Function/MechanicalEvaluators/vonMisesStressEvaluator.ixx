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
  real vonMisesStressEvaluator<N, StressEvaluatorType>::apply(
      const auto& values) const {
    using namespace tfel::math;
    const auto sig = map<stensor<N, real>>(values.data());
    return sigmaeq(sig);
  }

  namespace internals {

    template <unsigned short N>
    template <typename StressEvaluatorType>
    constexpr auto vmis_modifier<N>::operator()(StressEvaluatorType&& e) const
        requires(EvaluatorConcept<std::decay_t<StressEvaluatorType>>) {
      return vonMisesStressEvaluator<N, std::decay_t<StressEvaluatorType>>(
          std::forward<StressEvaluatorType>(e));
    }  // end of vmis

  }  // namespace internals

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_VONMISESSTRESSEVALUATOR_IXX */
