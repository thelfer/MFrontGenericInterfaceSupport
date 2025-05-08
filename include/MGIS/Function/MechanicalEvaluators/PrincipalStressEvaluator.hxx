/*!
 * \file   PrincipalStressEvaluator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   01/05/2025
 */

#ifndef MGIS_HAVE_TFEL
#error "TFEL is required to use mechanical evaluators"
#endif /* MGIS_HAVE_TFEL */

#ifndef LIB_MGIS_FUNCTION_PRINCIPALSTRESSEVALUATOR_HXX
#define LIB_MGIS_FUNCTION_PRINCIPALSTRESSEVALUATOR_HXX

#include <array>
#include <TFEL/Math/stensor.hxx>
#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/MechanicalEvaluators/StressEvaluatorBase.hxx"

namespace mgis::function {

  /*!
   * \brief an evaluator returning the von Mises stress from
   * an evaluator of a stress, assumed to be a symmetric tensor
   * \tparam N: space dimension
   * \tparam StressEvaluatorType: evaluator of the stress
   * \tparam esolver: eigen solver to be used
   */
  template <unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            tfel::math::stensor_common::EigenSolver esolver =
                tfel::math::stensor_common::TFELEIGENSOLVER>
  requires((N == 1) || (N == 2) || (N == 3)) struct PrincipalStressEvaluator
      : StressEvaluatorBase<
            PrincipalStressEvaluator<N, StressEvaluatorType, esolver>,
            N,
            StressEvaluatorType,
            true> {
    // constructors
    using StressEvaluatorBase<
        PrincipalStressEvaluator<N, StressEvaluatorType, esolver>,
        N,
        StressEvaluatorType,
        true>::StressEvaluatorBase;
    //! \return the number of components
    constexpr size_type getNumberOfComponents() const noexcept;
    /*!
     * \brief call operator
     * \param[in] values: stress values
     */
    std::array<real, 3u> apply(const auto&) const;
  };

  template <unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            tfel::math::stensor_common::EigenSolver esolver =
                tfel::math::stensor_common::TFELEIGENSOLVER>
  auto principal_stress(const StressEvaluatorType&) requires((N == 1) ||
                                                             (N == 2) ||
                                                             (N == 3));

}  // namespace mgis::function

#include "MGIS/Function/MechanicalEvaluators/PrincipalStressEvaluator.ixx"

#endif /* LIB_MGIS_FUNCTION_PRINCIPALSTRESSEVALUATOR_HXX */
