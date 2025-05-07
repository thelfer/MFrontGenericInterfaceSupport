/*!
 * \file   MGIS/Function/MechanicalEvaluators/StressEvaluatorBase.ixx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_STRESSEVALUATORBASE_IXX
#define LIB_MGIS_FUNCTION_STRESSEVALUATORBASE_IXX

#include "TFEL/Math/tensor.hxx"
#include "TFEL/Math/stensor.hxx"

namespace mgis::function {

  template <typename Child,
            unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  bool StressEvaluatorBase<Child, N, StressEvaluatorType, isSymmetric>::check(
      Context& ctx) const {
    using namespace tfel::math;
    if (!EvaluatorModifierBase<Child, StressEvaluatorType>::check(ctx)) {
      return false;
    }
    const auto nc = this->evaluator.getNumberOfComponents();
    if constexpr (isSymmetric) {
      if (nc != StensorDimeToSize<N>::value) {
        return ctx.registerErrorMessage(
            "incompatible number of components of the stress");
      }
    } else {
      if (nc != TensorDimeToSize<N>::value) {
        return ctx.registerErrorMessage(
            "incompatible number of components of the stress");
      }
    }
    return true;
  }  // end of check

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_STRESSEVALUATORBASE_IXX */
