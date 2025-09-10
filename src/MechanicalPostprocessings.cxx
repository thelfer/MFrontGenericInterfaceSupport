/*!
 * \file   MechanicalPostprocessings.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   01/09/2025
 */

#ifdef MGIS_HAVE_TFEL

#include <optional>
#include "MGIS/AbstractErrorHandler.hxx"
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/Algorithms.hxx"
#include "MGIS/Function/Tensors.hxx"
#include "MGIS/Function/Mechanics.hxx"

namespace mgis::function {

  namespace internals {

    template <EvaluatorConcept EvaluatorType>
    std::optional<Function<BasicLinearSpace>> evaluate(
        AbstractErrorHandler& ctx, const EvaluatorType& e) {
      const auto& s = getSpace(e);
      const auto nc = getNumberOfComponents(e);
      if (!Function<BasicLinearSpace>::checkPreconditions(ctx, s, nc)) {
        return {};
      }
      auto r = Function<BasicLinearSpace>{s, nc};
      if (!assign(ctx, r, e)) {
        return {};
      }
      return r;
    }  // end of evaluate

    template <EvaluatorModifierConcept ModifierType>
    std::optional<Function<BasicLinearSpace>>
    evaluateUnaryOperationOnSymmetricTensor(
        AbstractErrorHandler& ctx,
        const FunctionView<BasicLinearSpace, {}, false>& f,
        const ModifierType& op) {
      if (getNumberOfComponents(f) == 3) {
        return internals::evaluate(ctx, f | as_stensor<1> | op);
      } else if (getNumberOfComponents(f) == 4) {
        return internals::evaluate(ctx, f | as_stensor<2> | op);
      } else if (getNumberOfComponents(f) == 6) {
        return internals::evaluate(ctx, f | as_stensor<3> | op);
      }
      return ctx.registerErrorMessage("invalid number of components");
    }  // end of evaluateUnaryOperationOnSymmetricTensor

  }  // end of namespace internals

  std::optional<Function<BasicLinearSpace>> computeTrace(
      AbstractErrorHandler& ctx,
      const FunctionView<BasicLinearSpace, {}, false>& f) {
    return internals::evaluateUnaryOperationOnSymmetricTensor(ctx, f, trace);
  }  // end of computeTrace

  std::optional<Function<BasicLinearSpace>> computeMisesStress(
      AbstractErrorHandler& ctx,
      const FunctionView<BasicLinearSpace, {}, false>& f) {
    return internals::evaluateUnaryOperationOnSymmetricTensor(ctx, f, vmis);
  }  // end of computeMisesStress

  std::optional<Function<BasicLinearSpace>> computeHydrostraticPressure(
      AbstractErrorHandler& ctx,
      const FunctionView<BasicLinearSpace, {}, false>& f) {
    return internals::evaluateUnaryOperationOnSymmetricTensor(
        ctx, f, hydrostatic_stress);
  }  // end of computeHydrostraticPressure

}  // end of namespace mgis::function

#endif /* MGIS_HAVE_TFEL */