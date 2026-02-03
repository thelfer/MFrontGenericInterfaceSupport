/*!
 * \file   MGIS/Function/TFEL/TensorOperations.hxx
 * \brief
 * \author Thomas Helfer
 * \date   15/11/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef MGIS_HAVE_TFEL
#error "TFEL is required to use tensor evaluators"
#endif /* MGIS_HAVE_TFEL */

#ifndef LIB_MGIS_FUNCTION_TFEL_TENSOROPERATIONS_HXX
#define LIB_MGIS_FUNCTION_TFEL_TENSOROPERATIONS_HXX

#include <span>
#include <optional>
#include "MGIS/Config.hxx"
#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/Function.hxx"

namespace mgis::function {

  /*!
   * \brief a simple alias to shorten the declaration of the functions
   */
  using BasicImmutableFunctionView =
      FunctionView<BasicLinearSpace, FunctionDataLayoutDescription{}, false>;

  /*!
   * \brief type of tensor that are supported
   */
  enum struct TensorType { STENSOR, TENSOR, ST2TOST2 };
  /*!
   * \brief supported rotation operations
   *
   * - If the rotation operation is `FORWARD`, the given rotation matrix is used
   *   directly.
   * - If the rotation operation is `BACKWARD`, the transposed of the rotation
   *   matrix is used.
   */
  enum struct RotationOperation { FORWARD, BACKWARD };

  MGIS_EXPORT std::optional<Function<BasicLinearSpace>> computeRotatedTensor(
      Context&,
      const BasicImmutableFunctionView&,
      const BasicImmutableFunctionView&,
      const TensorType = TensorType::STENSOR,
      const RotationOperation = RotationOperation::FORWARD) noexcept;

  MGIS_EXPORT std::optional<Function<BasicLinearSpace>> computeRotatedTensor(
      Context&,
      const BasicImmutableFunctionView&,
      const std::span<real>&,
      const TensorType = TensorType::STENSOR,
      const RotationOperation = RotationOperation::FORWARD) noexcept;

  MGIS_EXPORT std::optional<Function<BasicLinearSpace>>
  computeRotatedSymmetricTensor(
      Context&,
      const BasicImmutableFunctionView&,
      const BasicImmutableFunctionView&,
      const RotationOperation = RotationOperation::FORWARD) noexcept;

  MGIS_EXPORT std::optional<Function<BasicLinearSpace>>
  computeRotatedSymmetricTensor(
      Context&,
      const BasicImmutableFunctionView&,
      const std::span<real>&,
      const RotationOperation = RotationOperation::FORWARD) noexcept;

  MGIS_EXPORT std::optional<Function<BasicLinearSpace>>
  computeRotatedUnsymmetricTensor(
      Context&,
      const BasicImmutableFunctionView&,
      const BasicImmutableFunctionView&,
      const RotationOperation = RotationOperation::FORWARD) noexcept;

  MGIS_EXPORT std::optional<Function<BasicLinearSpace>>
  computeRotatedUnsymmetricTensor(
      Context&,
      const BasicImmutableFunctionView&,
      const std::span<real>&,
      const RotationOperation = RotationOperation::FORWARD) noexcept;

  MGIS_EXPORT std::optional<Function<BasicLinearSpace>> computeRotatedST2toST2(
      Context&,
      const BasicImmutableFunctionView&,
      const BasicImmutableFunctionView&,
      const RotationOperation = RotationOperation::FORWARD) noexcept;

  MGIS_EXPORT std::optional<Function<BasicLinearSpace>> computeRotatedST2toST2(
      Context&,
      const BasicImmutableFunctionView&,
      const std::span<real>&,
      const RotationOperation = RotationOperation::FORWARD) noexcept;

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TFEL_TENSOROPERATIONS_HXX */
