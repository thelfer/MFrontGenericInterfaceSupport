/*!
 * \file   MGIS/Function/TFEL/TensorOperations.hxx
 * \brief  File which declares functions rotating the tensorial values of a
 * function into another basis and returning the result as a new function.
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
#include "MGIS/Function/TFEL/General.hxx"

namespace mgis::function {

  /*!
   * \brief immutable view on a function defined on a `BasicLinearSpace` whose
   * number of components is only known at run time
   *
   * This alias shortens the declaration of the functions of this file. The
   * `false` template argument states that the viewed values cannot be modified.
   */
  using BasicImmutableFunctionView =
      FunctionView<BasicLinearSpace, FunctionDataLayoutDescription{}, false>;

  /*!
   * \brief type of tensor that are supported
   *
   * - `STENSOR` designates a symmetric tensor,
   * - `TENSOR` designates an unsymmetric tensor,
   * - `ST2TOST2` designates a fourth order tensor mapping symmetric tensors to
   *   symmetric tensors.
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

  /*!
   * \brief rotate the tensorial values of a function and return the result as a
   * new function
   *
   * The type of the tensorial objects handled is given by the `type` argument
   * and the rotation is performed with the rotation matrix given by the `r`
   * view. The number of components of `t` must be compatible with `type`:
   * 3, 4 or 6 components for a symmetric tensor, 3, 5 or 9 components for an
   * unsymmetric tensor and 9, 16 or 36 components for a fourth order tensor
   * `st2tost2`.
   *
   * \param[in] ctx: execution context used to report the errors
   * \param[in] t: view on the function whose values are rotated
   * \param[in] r: view on the function whose values are the rotation matrices,
   * which must have 9 components
   * \param[in] type: type of the tensorial objects handled by the function
   * \param[in] ro: rotation operation to be performed
   * \return the function holding the rotated values or an empty optional if the
   * inputs are inconsistent, in which case an error is reported through the
   * error handler
   */
  MGIS_EXPORT std::optional<Function<BasicLinearSpace>> computeRotatedTensor(
      Context&,
      const BasicImmutableFunctionView&,
      const BasicImmutableFunctionView&,
      const TensorType = TensorType::STENSOR,
      const RotationOperation = RotationOperation::FORWARD) MGIS_NOEXCEPT;

  /*!
   * \brief rotate the tensorial values of a function with a given rotation
   * matrix and return the result as a new function
   *
   * This is the same operation as the previous overload, except that the
   * rotation matrix is given directly as a span of 9 values instead of being
   * taken from a function.
   *
   * \param[in] ctx: execution context used to report the errors
   * \param[in] t: view on the function whose values are rotated
   * \param[in] r: values of the rotation matrix, which must have 9 components
   * \param[in] type: type of the tensorial objects handled by the function
   * \param[in] ro: rotation operation to be performed
   * \return the function holding the rotated values or an empty optional if the
   * inputs are inconsistent, in which case an error is reported through the
   * error handler
   */
  MGIS_EXPORT std::optional<Function<BasicLinearSpace>> computeRotatedTensor(
      Context&,
      const BasicImmutableFunctionView&,
      const std::span<real>&,
      const TensorType = TensorType::STENSOR,
      const RotationOperation = RotationOperation::FORWARD) MGIS_NOEXCEPT;

  /*!
   * \brief rotate the values of a function seen as symmetric tensors and return
   * the result as a new function
   *
   * The number of components of `s` must be 3, 4 or 6, corresponding to
   * symmetric tensors of dimension 1, 2 or 3. The rotation matrix is taken from
   * the `r` view, which must have 9 components.
   *
   * \param[in] ctx: execution context used to report the errors
   * \param[in] s: view on the function whose symmetric tensorial values are
   * rotated
   * \param[in] r: view on the function whose values are the rotation matrices
   * \param[in] ro: rotation operation to be performed
   * \return the function holding the rotated values or an empty optional if the
   * inputs are inconsistent, in which case an error is reported through the
   * error handler
   */
  MGIS_EXPORT std::optional<Function<BasicLinearSpace>>
  computeRotatedSymmetricTensor(
      Context&,
      const BasicImmutableFunctionView&,
      const BasicImmutableFunctionView&,
      const RotationOperation = RotationOperation::FORWARD) MGIS_NOEXCEPT;

  /*!
   * \brief rotate the values of a function seen as symmetric tensors with a
   * given rotation matrix and return the result as a new function
   *
   * This is the same operation as the previous overload, except that the
   * rotation matrix is given directly as a span of 9 values.
   *
   * \param[in] ctx: execution context used to report the errors
   * \param[in] s: view on the function whose symmetric tensorial values are
   * rotated
   * \param[in] r: values of the rotation matrix, which must have 9 components
   * \param[in] ro: rotation operation to be performed
   * \return the function holding the rotated values or an empty optional if the
   * inputs are inconsistent, in which case an error is reported through the
   * error handler
   */
  MGIS_EXPORT std::optional<Function<BasicLinearSpace>>
  computeRotatedSymmetricTensor(
      Context&,
      const BasicImmutableFunctionView&,
      const std::span<real>&,
      const RotationOperation = RotationOperation::FORWARD) MGIS_NOEXCEPT;

  /*!
   * \brief rotate the values of a function seen as unsymmetric tensors and
   * return the result as a new function
   *
   * The number of components of `s` must be 3, 5 or 9, corresponding to
   * unsymmetric tensors of dimension 1, 2 or 3. The rotation matrix is taken
   * from the `r` view, which must have 9 components.
   *
   * \param[in] ctx: execution context used to report the errors
   * \param[in] s: view on the function whose tensorial values are rotated
   * \param[in] r: view on the function whose values are the rotation matrices
   * \param[in] ro: rotation operation to be performed
   * \return the function holding the rotated values or an empty optional if the
   * inputs are inconsistent, in which case an error is reported through the
   * error handler
   */
  MGIS_EXPORT std::optional<Function<BasicLinearSpace>>
  computeRotatedUnsymmetricTensor(
      Context&,
      const BasicImmutableFunctionView&,
      const BasicImmutableFunctionView&,
      const RotationOperation = RotationOperation::FORWARD) MGIS_NOEXCEPT;

  /*!
   * \brief rotate the values of a function seen as unsymmetric tensors with a
   * given rotation matrix and return the result as a new function
   *
   * This is the same operation as the previous overload, except that the
   * rotation matrix is given directly as a span of 9 values.
   *
   * \param[in] ctx: execution context used to report the errors
   * \param[in] s: view on the function whose tensorial values are rotated
   * \param[in] r: values of the rotation matrix, which must have 9 components
   * \param[in] ro: rotation operation to be performed
   * \return the function holding the rotated values or an empty optional if the
   * inputs are inconsistent, in which case an error is reported through the
   * error handler
   */
  MGIS_EXPORT std::optional<Function<BasicLinearSpace>>
  computeRotatedUnsymmetricTensor(
      Context&,
      const BasicImmutableFunctionView&,
      const std::span<real>&,
      const RotationOperation = RotationOperation::FORWARD) MGIS_NOEXCEPT;

  /*!
   * \brief rotate the values of a function seen as fourth order tensors
   * mapping symmetric tensors to symmetric tensors and return the result as a
   * new function
   *
   * The number of components of `s` must be 9, 16 or 36, corresponding to
   * fourth order tensors of dimension 1, 2 or 3. The rotation matrix is taken
   * from the `r` view, which must have 9 components.
   *
   * \param[in] ctx: execution context used to report the errors
   * \param[in] s: view on the function whose fourth order tensorial values are
   * rotated
   * \param[in] r: view on the function whose values are the rotation matrices
   * \param[in] ro: rotation operation to be performed
   * \return the function holding the rotated values or an empty optional if the
   * inputs are inconsistent, in which case an error is reported through the
   * error handler
   */
  MGIS_EXPORT std::optional<Function<BasicLinearSpace>> computeRotatedST2toST2(
      Context&,
      const BasicImmutableFunctionView&,
      const BasicImmutableFunctionView&,
      const RotationOperation = RotationOperation::FORWARD) MGIS_NOEXCEPT;

  /*!
   * \brief rotate the values of a function seen as fourth order tensors
   * mapping symmetric tensors to symmetric tensors with a given rotation matrix
   * and return the result as a new function
   *
   * This is the same operation as the previous overload, except that the
   * rotation matrix is given directly as a span of 9 values.
   *
   * \param[in] ctx: execution context used to report the errors
   * \param[in] s: view on the function whose fourth order tensorial values are
   * rotated
   * \param[in] r: values of the rotation matrix, which must have 9 components
   * \param[in] ro: rotation operation to be performed
   * \return the function holding the rotated values or an empty optional if the
   * inputs are inconsistent, in which case an error is reported through the
   * error handler
   */
  MGIS_EXPORT std::optional<Function<BasicLinearSpace>> computeRotatedST2toST2(
      Context&,
      const BasicImmutableFunctionView&,
      const std::span<real>&,
      const RotationOperation = RotationOperation::FORWARD) MGIS_NOEXCEPT;

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TFEL_TENSOROPERATIONS_HXX */
