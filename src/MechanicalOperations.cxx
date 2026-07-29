/*!
 * \file   src/MechanicalOperations.cxx
 * \brief
 * \author Thomas Helfer
 * \date   01/09/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <optional>
#include "MGIS/Context.hxx"
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/Algorithms.hxx"
#include "MGIS/Function/TFEL/Tensors.hxx"
#include "MGIS/Function/TFEL/Mechanics.hxx"
#include "MGIS/Function/TFEL/MechanicalOperations.hxx"

namespace mgis::function::internals {

  template <EvaluatorConcept EvaluatorType>
  static std::optional<Function<BasicLinearSpace>> evaluate(
      Context& ctx, const EvaluatorType& e) noexcept {
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

  static std::optional<Function<BasicLinearSpace>>
  evaluateUnaryOperationOnSymmetricTensor(
      Context& ctx,
      const BasicImmutableFunctionView& f,
      const EvaluatorModifierConcept auto& op) noexcept {
    if (getNumberOfComponents(f) == 3) {
      return internals::evaluate(ctx, f | as_stensor<1> | op);
    } else if (getNumberOfComponents(f) == 4) {
      return internals::evaluate(ctx, f | as_stensor<2> | op);
    } else if (getNumberOfComponents(f) == 6) {
      return internals::evaluate(ctx, f | as_stensor<3> | op);
    }
    return ctx.registerErrorMessage("invalid number of components");
  }  // end of evaluateUnaryOperationOnSymmetricTensor

  static std::optional<Function<BasicLinearSpace>>
  evaluateUnaryOperationOnCauchyStressComputedFromPK1(
      Context& ctx,
      const BasicImmutableFunctionView& pk1,
      const BasicImmutableFunctionView& F,
      const EvaluatorModifierConcept auto& op) noexcept {
    if (getNumberOfComponents(pk1) == 3) {
      return internals::evaluate(
          ctx, pk1 | as_tensor<1> | from_pk1_to_cauchy(F | as_tensor<1>) | op);
    } else if (getNumberOfComponents(pk1) == 5) {
      return internals::evaluate(
          ctx, pk1 | as_tensor<2> | from_pk1_to_cauchy(F | as_tensor<2>) | op);
    } else if (getNumberOfComponents(pk1) == 9) {
      return internals::evaluate(
          ctx, pk1 | as_tensor<3> | from_pk1_to_cauchy(F | as_tensor<3>) | op);
    }
    return ctx.registerErrorMessage("invalid number of components");
  }  // end of evaluateUnaryOperationOnUnsymmetricTensor

}  // namespace mgis::function::internals

namespace mgis::function {

  std::optional<Function<BasicLinearSpace>> computeTrace(
      Context& ctx, const BasicImmutableFunctionView& f) noexcept {
    return internals::evaluateUnaryOperationOnSymmetricTensor(ctx, f, trace);
  }  // end of computeTrace

  std::optional<Function<BasicLinearSpace>> computeVonMisesStress(
      Context& ctx, const BasicImmutableFunctionView& f) noexcept {
    return internals::evaluateUnaryOperationOnSymmetricTensor(ctx, f, vmis);
  }  // end of computeMisesStress

  std::optional<Function<BasicLinearSpace>> computeHydrostraticPressure(
      Context& ctx, const BasicImmutableFunctionView& f) noexcept {
    return internals::evaluateUnaryOperationOnSymmetricTensor(
        ctx, f, hydrostatic_stress);
  }  // end of computeHydrostraticPressure

  std::optional<Function<BasicLinearSpace>> computeEigenValues(
      Context& ctx, const BasicImmutableFunctionView& f) noexcept {
    return internals::evaluateUnaryOperationOnSymmetricTensor(ctx, f,
                                                              eigen_values<>);
  }  // end of computeHydrostraticPressure

  std::optional<Function<BasicLinearSpace>>
  computeCauchyStressFromFirstPiolaKirchhoffStress(
      Context& ctx,
      const BasicImmutableFunctionView& pk1,
      const BasicImmutableFunctionView& F) noexcept {
    if (getNumberOfComponents(pk1) == 3) {
      return internals::evaluate(
          ctx, pk1 | as_tensor<1> | from_pk1_to_cauchy(F | as_tensor<1>));
    } else if (getNumberOfComponents(pk1) == 5) {
      return internals::evaluate(
          ctx, pk1 | as_tensor<2> | from_pk1_to_cauchy(F | as_tensor<2>));
    } else if (getNumberOfComponents(pk1) == 9) {
      return internals::evaluate(
          ctx, pk1 | as_tensor<3> | from_pk1_to_cauchy(F | as_tensor<3>));
    }
    return ctx.registerErrorMessage("invalid number of components");
  }  // end of computeCauchyStressFromFirstPiolaKirchhoffStress

  std::optional<Function<BasicLinearSpace>>
  computeCauchyStressHydrostraticPressureFromFirstPiolaKirchhoffStress(
      Context& ctx,
      const BasicImmutableFunctionView& pk1,
      const BasicImmutableFunctionView& F) noexcept {
    return internals::evaluateUnaryOperationOnCauchyStressComputedFromPK1(
        ctx, pk1, F, vmis);
  }

  std::optional<Function<BasicLinearSpace>>
  computeCauchyStressVonMisesStressFromFirstPiolaKirchhoffStress(
      Context& ctx,
      const BasicImmutableFunctionView& pk1,
      const BasicImmutableFunctionView& F) noexcept {
    return internals::evaluateUnaryOperationOnCauchyStressComputedFromPK1(
        ctx, pk1, F, vmis);
  }  // end of computeCauchyStressVonMisesStressFromFirstPiolaKirchhoffStress

  std::optional<Function<BasicLinearSpace>>
  computeCauchyStressEigenValuesFromFirstPiolaKirchhoffStress(
      Context& ctx,
      const BasicImmutableFunctionView& pk1,
      const BasicImmutableFunctionView& F) noexcept {
    return internals::evaluateUnaryOperationOnCauchyStressComputedFromPK1(
        ctx, pk1, F, eigen_values<>);
  }  // end of computeStressCauchyEigenValuesFromFirstPiolaKirchhoffStress

}  // end of namespace mgis::function
