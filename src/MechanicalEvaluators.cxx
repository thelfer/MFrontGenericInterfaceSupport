/*!
 * \file   MechanicalEvaluators.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#include "MGIS/QuadratureFunction/Evaluators.hxx"
#include "MGIS/QuadratureFunction/MechanicalEvaluators.hxx"

namespace mgis::quadrature_function {

#ifdef MGIS_HAVE_TFEL

  static_assert(
      EvaluatorConcept<CauchyStressFromFirstPiolaKirchhoffStressEvaluator<
          3,
          FixedSizedEvaluator<9>,
          FixedSizedEvaluator<9>>>);
  static_assert(
      EvaluatorConcept<vonMisesStressEvaluator<3, FixedSizedEvaluator<6>>>);
  static_assert(
      EvaluatorConcept<PrincipalStressEvaluator<3, FixedSizedEvaluator<6>>>);

#endif MGIS_HAVE_TFEL

}  // end of namespace mgis::quadrature_function