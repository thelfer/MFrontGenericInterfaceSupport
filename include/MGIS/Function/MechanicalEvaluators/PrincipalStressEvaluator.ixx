/*!
 * \file   PrincipalStressEvaluator.ixx
 * \brief
 * \author Thomas Helfer
 * \date   01/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_PRINCIPALSTRESSEVALUATOR_IXX
#define LIB_MGIS_FUNCTION_PRINCIPALSTRESSEVALUATOR_IXX

#include "TFEL/Math/stensor.hxx"
#include "TFEL/Math/Array/View.hxx"

namespace mgis::function {

  template <unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            tfel::math::stensor_common::EigenSolver esolver>
  constexpr size_type PrincipalStressEvaluator<N,
                                               StressEvaluatorType,
                                               esolver>::getNumberOfComponents()
      const noexcept {
    return 3;
  }  // end of getNumberOfComponents

  template <unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            tfel::math::stensor_common::EigenSolver esolver>
  std::array<real, 3u>
  PrincipalStressEvaluator<N, StressEvaluatorType, esolver>::apply(
      const auto& values) const {
    using namespace tfel::math;
    const auto sig = stensor<N, real>(values.data());
    const auto vp = sig.template computeEigenValues<esolver>();
    return {vp[0], vp[1], vp[2]};
  }  // end of operator()

  namespace internals {

    template <unsigned short N, tfel::math::stensor_common::EigenSolver esolver>
    template <typename StressEvaluatorType>
    constexpr auto principal_stress_modifier<N, esolver>::operator()(
        StressEvaluatorType&& e) const
        requires(EvaluatorConcept<std::decay_t<StressEvaluatorType>>) {
      return PrincipalStressEvaluator<N, std::decay_t<StressEvaluatorType>,
                                      esolver>(
          std::forward<StressEvaluatorType>(e));
    }  // end of principal_stress

  }  // namespace internals

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_PRINCIPALSTRESSEVALUATOR_IXX */
