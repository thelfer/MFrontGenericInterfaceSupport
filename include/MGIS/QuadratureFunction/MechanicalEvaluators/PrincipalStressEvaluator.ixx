/*!
 * \file   PrincipalStressEvaluator.ixx
 * \brief
 * \author Thomas Helfer
 * \date   01/05/2025
 */

#ifndef LIB_MGIS_QUADRATUREFUNCTION_PRINCIPALSTRESSEVALUATOR_IXX
#define LIB_MGIS_QUADRATUREFUNCTION_PRINCIPALSTRESSEVALUATOR_IXX

#include "TFEL/Math/stensor.hxx"
#include "TFEL/Math/Array/View.hxx"

namespace mgis::quadrature_function {

  template <unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            tfel::math::stensor_common::EigenSolver esolver>
  PrincipalStressEvaluator<N, StressEvaluatorType, esolver>::
      PrincipalStressEvaluator(const StressEvaluatorType& e, const real eps)
      : StressEvaluatorBase<N, StressEvaluatorType, true>(e), seps(eps) {}

  template <unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            tfel::math::stensor_common::EigenSolver esolver>
  PrincipalStressEvaluator<N, StressEvaluatorType, esolver>::
      PrincipalStressEvaluator(const PrincipalStressEvaluator&) = default;

  template <unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            tfel::math::stensor_common::EigenSolver esolver>
  PrincipalStressEvaluator<N, StressEvaluatorType, esolver>::
      PrincipalStressEvaluator(PrincipalStressEvaluator&&) = default;

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
  PrincipalStressEvaluator<N, StressEvaluatorType, esolver>::operator()(
      const size_type i) const {
    using namespace tfel::math;
    const auto sig = stensor<N, real>(this->stress_evaluator(i).data());
    const auto vp = sig.template computeEigenValues<esolver>(this->seps);
    return {vp[0], vp[1], vp[2]};
  }  // end of operator()

}  // namespace mgis::quadrature_function

#endif /* LIB_MGIS_QUADRATUREFUNCTION_PRINCIPALSTRESSEVALUATOR_IXX */
