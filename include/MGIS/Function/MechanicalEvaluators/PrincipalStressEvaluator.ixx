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

  template <unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            tfel::math::stensor_common::EigenSolver esolver =
                tfel::math::stensor_common::TFELEIGENSOLVER>
  auto principal_stress(const StressEvaluatorType& e) requires((N == 1) ||
                                                               (N == 2) ||
                                                               (N == 3)) {
    return PrincipalStressEvaluator<N, StressEvaluatorType, esolver>(e);

  }  // end of principal_stress

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_PRINCIPALSTRESSEVALUATOR_IXX */
