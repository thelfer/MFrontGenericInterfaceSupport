/*!
 * \file   src/TensorOperations.cxx
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
#include "MGIS/Function/TFEL/TensorOperations.hxx"

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

  template <typename TensorFunctionView>
  std::optional<Function<BasicLinearSpace>> computeRotatedTensor(
      Context& ctx,
      const TensorFunctionView& t,
      const BasicImmutableFunctionView& r,
      const RotationOperation ro) noexcept {
    const auto rm = r | as_tmatrix<3, 3>;
    if (ro == RotationOperation::FORWARD) {
      return internals::evaluate(ctx, t | rotate(rm));
    } else {
      return internals::evaluate(ctx, t | rotate_backwards(rm));
    }
    return ctx.registerErrorMessage("invalid number of components");
  }  // end of computeRotatedSymmetricTensor

  template <typename TensorFunctionView>
  std::optional<Function<BasicLinearSpace>> computeRotatedTensor(
      Context& ctx,
      const TensorFunctionView& t,
      const std::span<real>& r,
      const RotationOperation ro) noexcept {
    const auto rm = tfel::math::tmatrix<3u, 3u, real>{r.data()};
    if (ro == RotationOperation::FORWARD) {
      return internals::evaluate(ctx, t | rotate(rm));
    } else {
      return internals::evaluate(ctx, t | rotate_backwards(rm));
    }
    return ctx.registerErrorMessage("invalid number of components");
  }  // end of computeRotatedSymmetricTensor

} // end of namespace mgis::function::internals

namespace mgis::function {

  std::optional<Function<BasicLinearSpace>> computeRotatedTensor(
      Context& ctx,
      const BasicImmutableFunctionView& t,
      const BasicImmutableFunctionView& r,
      const TensorType type,
      const RotationOperation ro) noexcept{
    if (type == TensorType::STENSOR) {
      return computeRotatedSymmetricTensor(ctx, t, r, ro);
    }
    return ctx.registerErrorMessage(
        "computeRotatedTensor: unsupported tensor type");
  }

  std::optional<Function<BasicLinearSpace>> computeRotatedTensor(
      Context& ctx,
      const BasicImmutableFunctionView& t,
      const std::span<real>& r,
      const TensorType type,
      const RotationOperation ro) noexcept {
    if (type == TensorType::STENSOR) {
      return computeRotatedSymmetricTensor(ctx, t, r, ro);
    } else if (type == TensorType::TENSOR) {
      return computeRotatedUnsymmetricTensor(ctx, t, r, ro);
    }
    return ctx.registerErrorMessage(
        "computeRotatedTensor: unsupported tensor type");
  } // end of computeRotatedTensor

  std::optional<Function<BasicLinearSpace>> computeRotatedSymmetricTensor(
      Context& ctx,
      const BasicImmutableFunctionView& s,
      const BasicImmutableFunctionView& r,
      const RotationOperation ro) noexcept {
    if (getNumberOfComponents(s) == 3) {
      return internals::evaluate(ctx, s);
    } else if (getNumberOfComponents(s) == 4) {
      return internals::computeRotatedTensor(ctx, s | as_stensor<2>, r, ro);
    } else if (getNumberOfComponents(s) == 6) {
      return internals::computeRotatedTensor(ctx, s | as_stensor<3>, r, ro);
    }
    return ctx.registerErrorMessage("invalid number of components");
  }  // end of computeRotatedSymmetricTensor

  std::optional<Function<BasicLinearSpace>> computeRotatedSymmetricTensor(
      Context& ctx,
      const BasicImmutableFunctionView& s,
      const std::span<real>& r,
      const RotationOperation ro) noexcept {
    if (r.size() != 9u) {
      return ctx.registerErrorMessage(
          "computeRotatedSymmetricTensor: invalid number of values for the "
          "rotation matrix");
    }
    if (getNumberOfComponents(s) == 3) {
      return internals::evaluate(ctx, s);
    } else if (getNumberOfComponents(s) == 4) {
      return internals::computeRotatedTensor(ctx, s | as_stensor<2>, r, ro);
    } else if (getNumberOfComponents(s) == 6) {
      return internals::computeRotatedTensor(ctx, s | as_stensor<3>, r, ro);
    }
    return ctx.registerErrorMessage("invalid number of components");
  }  // end of computeRotatedSymmetricTensor

  std::optional<Function<BasicLinearSpace>> computeRotatedUnsymmetricTensor(
      Context& ctx,
      const BasicImmutableFunctionView& s,
      const BasicImmutableFunctionView& r,
      const RotationOperation ro) noexcept {
    if (getNumberOfComponents(s) == 3) {
      return internals::evaluate(ctx, s);
    } else if (getNumberOfComponents(s) == 5) {
      return internals::computeRotatedTensor(ctx, s | as_tensor<2>, r, ro);
    } else if (getNumberOfComponents(s) == 9) {
      return internals::computeRotatedTensor(ctx, s | as_tensor<3>, r, ro);
    }
    return ctx.registerErrorMessage("invalid number of components");
  }  // end of computeRotatedUnsymmetricTensor

  std::optional<Function<BasicLinearSpace>> computeRotatedUnsymmetricTensor(
      Context& ctx,
      const BasicImmutableFunctionView& s,
      const std::span<real>& r,
      const RotationOperation ro) noexcept {
    if (r.size() != 9u) {
      return ctx.registerErrorMessage(
          "computeRotatedUnsymmetricTensor: invalid number of values for the "
          "rotation matrix");
    }
    if (getNumberOfComponents(s) == 3) {
      return internals::evaluate(ctx, s);
    } else if (getNumberOfComponents(s) == 5) {
      return internals::computeRotatedTensor(ctx, s | as_tensor<2>, r, ro);
    } else if (getNumberOfComponents(s) == 9) {
      return internals::computeRotatedTensor(ctx, s | as_tensor<3>, r, ro);
    }
    return ctx.registerErrorMessage("invalid number of components");
  }  // end of computeRotatedUnsymmetricTensor

}  // end of namespace mgis::function
