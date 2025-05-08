/*!
 * \file   HydrostaticStressEvaluator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   01/05/2025
 */

#ifndef MGIS_HAVE_TFEL
#error "TFEL is required to use mechanical evaluators"
#endif /* MGIS_HAVE_TFEL */

#ifndef LIB_MGIS_FUNCTION_HYDROSTATICSTRESSSTRESSEVALUATOR_HXX
#define LIB_MGIS_FUNCTION_HYDROSTATICSTRESSSTRESSEVALUATOR_HXX

#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/MechanicalEvaluators/StressEvaluatorBase.hxx"

namespace mgis::function {

  /*!
   * \brief an evaluator returning the von Mises stress from
   * an evaluator of a stress, assumed to be a symmetric tensor
   * \tparam N: space dimension
   * \tparam StressEvaluatorType: evaluator of the stress
   */
  template <unsigned short N, EvaluatorConcept StressEvaluatorType>
  requires((N == 1) || (N == 2) || (N == 3)) struct HydrostaticStressEvaluator
      : StressEvaluatorBase<HydrostaticStressEvaluator<N, StressEvaluatorType>,
                            N,
                            StressEvaluatorType,
                            true> {
    // inheriting constructors
    using StressEvaluatorBase<HydrostaticStressEvaluator,
                              N,
                              StressEvaluatorType,
                              true>::StressEvaluatorBase;
    //! \return the number of components
    constexpr size_type getNumberOfComponents() const noexcept;
    /*!
     * \brief call operator
     * \param[in] values: values of the stress
     */
    real apply(const auto&) const;
  };

  template <unsigned short N, typename StressEvaluatorType>
  auto hydrostatic_stress(StressEvaluatorType&&) requires(
      EvaluatorConcept<std::decay_t<StressEvaluatorType>> &&
      ((N == 1) || (N == 2) || (N == 3)));

}  // namespace mgis::function

#include "MGIS/Function/MechanicalEvaluators/HydrostaticStressEvaluator.ixx"

#endif /* LIB_MGIS_FUNCTION_HYDROSTATICSTRESSSTRESSEVALUATOR_HXX */
