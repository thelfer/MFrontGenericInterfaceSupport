/*!
 * \file   MGIS/Function/MechanicalEvaluators/StressEvaluatorBase.hxx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#ifndef MGIS_HAVE_TFEL
#error "TFEL is required to use mechanical evaluators"
#endif /* MGIS_HAVE_TFEL */

#ifndef LIB_MGIS_FUNCTION_STRESSEVALUATORBASE_HXX
#define LIB_MGIS_FUNCTION_STRESSEVALUATORBASE_HXX

#include "MGIS/Function/Space.hxx"
#include "MGIS/Function/Evaluators.hxx"
#include "MGIS/Function/EvaluatorModifierBase.hxx"

namespace mgis::function {

  /*!
   * \brief a base class for evaluators modifying a stress tensor
   * \tparam Child: child class
   * \tparam Space: discretization space
   * \tparam N: space dimension
   * \tparam StressEvaluatorType: evaluator of the stress
   * \tparam symmetric: boolean stating of the stress tensor is symmetric
   */
  template <typename Child,
            unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  struct StressEvaluatorBase
      : EvaluatorModifierBase<Child, StressEvaluatorType> {
    // inheriting constructors
    using EvaluatorModifierBase<Child,
                                StressEvaluatorType>::EvaluatorModifierBase;
    //! \brief perform consistency checks
    void check() const;
  };

}  // namespace mgis::function

#include "MGIS/Function/MechanicalEvaluators/StressEvaluatorBase.ixx"

#endif /* LIB_MGIS_FUNCTION_STRESSEVALUATORBASE_HXX */
