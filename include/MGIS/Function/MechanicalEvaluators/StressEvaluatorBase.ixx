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
  void StressEvaluatorBase<Child, N, StressEvaluatorType, isSymmetric>::check()
      const {
    using namespace tfel::math;
    EvaluatorModifierBase<Child, StressEvaluatorType>::check();
    const auto nc = this->evaluator.getNumberOfComponents();
    if constexpr (isSymmetric) {
      raise_if(nc != StensorDimeToSize<N>::value,
               "StressEvaluatorBase::check: "
               "incompatible number of components of the stress");
    } else {
      raise_if(nc != TensorDimeToSize<N>::value,
               "StressEvaluatorBase::check: "
               "incompatible number of components of the stress");
    }
  }  // end of check

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_STRESSEVALUATORBASE_IXX */
