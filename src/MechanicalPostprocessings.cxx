/*!
 * \file   MechanicalPostprocessings.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   01/09/2025
 */

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

    template <bool is_mutable>
    static std::optional<FunctionView<BasicLinearSpace, {}, is_mutable>>
    view_impl(
        AbstractErrorHandler& ctx,
        typename FunctionView<BasicLinearSpace, {}, is_mutable>::ExternalData
            values,
        const FunctionDataLayout<{}>& l) {
      using FunctionView = FunctionView<BasicLinearSpace, {}, is_mutable>;
      const auto n = [&values, &l]() -> size_type {
        auto dv =
            std::div(static_cast<long long>(values.size()), l.getDataStride());
        return static_cast<size_type>(dv.quot);
      }();
      if (n == 0) {
        return ctx.registerErrorMessage("invalid values size");
      }
      const auto space = BasicLinearSpace{n};
      if (!FunctionView::checkPreconditions(ctx, space, values, l)) {
        return {};
      }
      return FunctionView{space, values, l};
    }  // end of view_impl

    static std::optional<FunctionView<BasicLinearSpace, {}, false>> view(
        AbstractErrorHandler& ctx,
        std::span<const real> values,
        const FunctionDataLayout<{}>& l) {
      return internals::view_impl<false>(ctx, values, l);
    }  // end of view

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
    evaluateUnaryOperationOnSymmetricTensor(AbstractErrorHandler& ctx,
                                            std::span<const real> values,
                                            const FunctionDataLayout<{}>& l,
                                            const ModifierType& op) {
      const auto s = view(ctx, values, l);
      if (isInvalid(s)) {
        return {};
      }
      if (l.getNumberOfComponents() == 3) {
        return internals::evaluate(ctx, *s | as_stensor<1> | op);
      } else if (l.getNumberOfComponents() == 4) {
        return internals::evaluate(ctx, *s | as_stensor<2> | op);
      } else if (l.getNumberOfComponents() == 6) {
        return internals::evaluate(ctx, *s | as_stensor<3> | op);
      }
      return ctx.registerErrorMessage("invalid number of components");
    }  // end of evaluateUnaryOperationOnSymmetricTensor

  }  // end of namespace internals

#ifdef MGIS_HAVE_TFEL

  std::optional<Function<BasicLinearSpace>> computeTrace(
      AbstractErrorHandler& ctx,
      std::span<const real> values,
      const FunctionDataLayout<{}>& l) {
    return internals::evaluateUnaryOperationOnSymmetricTensor(ctx, values, l,
                                                              trace);
  }  // end of computeTrace

  std::optional<Function<BasicLinearSpace>> computeMisesStress(
      AbstractErrorHandler& ctx,
      std::span<const real> values,
      const FunctionDataLayout<{}>& l) {
    return internals::evaluateUnaryOperationOnSymmetricTensor(ctx, values, l,
                                                              vmis);
  }  // end of computeMisesStress

  std::optional<Function<BasicLinearSpace>> computeHydrostraticPressure(
      AbstractErrorHandler& ctx,
      std::span<const real> values,
      const FunctionDataLayout<{}>& l) {
    return internals::evaluateUnaryOperationOnSymmetricTensor(
        ctx, values, l, hydrostatic_stress);
  }  // end of computeHydrostraticPressure

#endif /* MGIS_HAVE_TFEL */

}  // end of namespace mgis::function
