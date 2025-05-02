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
#include "MGIS/Function/Evaluators.hxx"
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
  struct PrincipalStressEvaluator
      : StressEvaluatorBase<N, StressEvaluatorType, true> {
    /*!
     * \brief default constructor
     * \param[in] e: stress evaluator
     * \param[in] eps: small value (relative to the stress magnitude)
     */
    PrincipalStressEvaluator(const StressEvaluatorType&, const real);
    //! \brief copy constructor
    PrincipalStressEvaluator(const PrincipalStressEvaluator&);
    //! \brief move constructor
    PrincipalStressEvaluator(PrincipalStressEvaluator&&);
    //! \return the number of components
    constexpr size_type getNumberOfComponents() const noexcept;
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    std::array<real, 3u> operator()(const size_type) const;

   private:
    /*!
     * \brief small value, relative to the expected stress values, required by
     * some eigensolvers.
     *
     * Typically 1e-15 * Young's modulus is a good estimate
     */
    const real seps;
  };

}  // namespace mgis::function

#include "MGIS/Function/MechanicalEvaluators/PrincipalStressEvaluator.ixx"

#endif /* LIB_MGIS_FUNCTION_PRINCIPALSTRESSEVALUATOR_HXX */
