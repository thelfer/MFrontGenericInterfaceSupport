/*!
 * \file   vonMisesStressEvaluator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   01/05/2025
 */

#ifndef MGIS_HAVE_TFEL
#error "TFEL is required to use mechanical evaluators"
#endif /* MGIS_HAVE_TFEL */

#ifndef LIB_MGIS_FUNCTION_VONMISESSTRESSEVALUATOR_HXX
#define LIB_MGIS_FUNCTION_VONMISESSTRESSEVALUATOR_HXX

#include "MGIS/Function/Evaluators.hxx"
#include "MGIS/Function/MechanicalEvaluators/StressEvaluatorBase.hxx"

namespace mgis::function {

  /*!
   * \brief an evaluator returning the von Mises stress from
   * an evaluator of a stress, assumed to be a symmetric tensor
   * \tparam N: space dimension
   * \tparam StressEvaluatorType: evaluator of the stress
   */
  template <unsigned short N, EvaluatorConcept StressEvaluatorType>
  requires((N == 1) || (N == 2) || (N == 3)) struct vonMisesStressEvaluator
      : StressEvaluatorBase<N, StressEvaluatorType, true> {
    // inheriting constructors
    using StressEvaluatorBase<N, StressEvaluatorType, true>::
        StressEvaluatorBase;
    //! \return the number of components
    constexpr size_type getNumberOfComponents() const noexcept;
    /*!
     * \brief call operator
     * \param[in] i: element index
     */
    real operator()(const size_type) const;
  };

  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))
  struct vmis_fn {
    template <typename StressEvaluatorType>
    constexpr auto operator()(StressEvaluatorType&&) const
        requires(EvaluatorConcept<std::decay_t<StressEvaluatorType>>);
  };

  template <unsigned short N>
  inline constexpr vmis_fn<N> vmis = {};

}  // namespace mgis::function

#include "MGIS/Function/MechanicalEvaluators/vonMisesStressEvaluator.ixx"

#endif /* LIB_MGIS_FUNCTION_VONMISESSTRESSEVALUATOR_HXX */
