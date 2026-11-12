/*!
 * \file   MGIS/Function/TFEL/MechanicalOperations.hxx
 * \brief  File which declares functions computing mechanical quantities (trace,
 * von Mises stress, hydrostatic pressure, eigen values, Cauchy stress) from the
 * tensorial values of a function and returning the result as a new function.
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

#ifndef LIB_MGIS_FUNCTION_TFEL_MECHANICALOPERATIONS_HXX
#define LIB_MGIS_FUNCTION_TFEL_MECHANICALOPERATIONS_HXX

#include <optional>
#include "MGIS/Config.hxx"
#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/Function.hxx"

namespace mgis::function {

  //! \brief immutable view on a function defined on a `BasicLinearSpace` whose
  //! number of components is only known at run time
  using BasicImmutableFunctionView =
      FunctionView<BasicLinearSpace, FunctionDataLayoutDescription{}, false>;

  /*!
   * \brief compute the trace of the symmetric tensorial values of a function
   *
   * The values of `f` are interpreted as symmetric tensors, which requires 3,
   * 4 or 6 components, corresponding to symmetric tensors of dimension 1, 2 or
   * 3; the result is a scalar valued function.
   *
   * \param[in] ctx: execution context used to report the errors
   * \param[in] f: view on the function whose symmetric tensorial values are
   * considered
   * \return the function holding the traces or an empty optional if the number
   * of components of `f` is not compatible with a symmetric tensor, in which
   * case an error is reported through the error handler
   */
  MGIS_EXPORT std::optional<Function<BasicLinearSpace>> computeTrace(
      Context&, const BasicImmutableFunctionView&) MGIS_NOEXCEPT;

  /*!
   * \brief compute the von Mises stress of the symmetric tensorial values of a
   * function
   *
   * The values of `f` are interpreted as symmetric tensors, which requires 3,
   * 4 or 6 components, corresponding to symmetric tensors of dimension 1, 2 or
   * 3; the result is a scalar valued function.
   *
   * \param[in] ctx: execution context used to report the errors
   * \param[in] f: view on the function whose symmetric tensorial values are
   * interpreted as stresses
   * \return the function holding the von Mises stresses or an empty optional if
   * the number of components of `f` is not compatible with a symmetric tensor,
   * in which case an error is reported through the error handler
   */
  MGIS_EXPORT std::optional<Function<BasicLinearSpace>> computeVonMisesStress(
      Context&, const BasicImmutableFunctionView&) MGIS_NOEXCEPT;

  /*!
   * \brief compute the hydrostatic pressure of the symmetric tensorial values
   * of a function
   *
   * The hydrostatic pressure is defined as `trace(s) / 3`. The values of `f`
   * are interpreted as symmetric tensors, which requires 3, 4 or 6 components,
   * corresponding to symmetric tensors of dimension 1, 2 or 3; the result is a
   * scalar valued function.
   *
   * \param[in] ctx: execution context used to report the errors
   * \param[in] f: view on the function whose symmetric tensorial values are
   * interpreted as stresses
   * \return the function holding the hydrostatic pressures or an empty optional
   * if the number of components of `f` is not compatible with a symmetric
   * tensor, in which case an error is reported through the error handler
   */
  MGIS_EXPORT std::optional<Function<BasicLinearSpace>>
  computeHydrostraticPressure(Context&,
                              const BasicImmutableFunctionView&) MGIS_NOEXCEPT;

  /*!
   * \brief compute the eigen values of the symmetric tensorial values of a
   * function
   *
   * The values of `f` are interpreted as symmetric tensors, which requires 3,
   * 4 or 6 components, corresponding to symmetric tensors of dimension 1, 2 or
   * 3; the result is a function whose values have as many components as the
   * dimension of the symmetric tensors.
   *
   * \param[in] ctx: execution context used to report the errors
   * \param[in] f: view on the function whose symmetric tensorial values are
   * considered
   * \return the function holding the eigen values or an empty optional if the
   * number of components of `f` is not compatible with a symmetric tensor, in
   * which case an error is reported through the error handler
   */
  MGIS_EXPORT std::optional<Function<BasicLinearSpace>> computeEigenValues(
      Context&, const BasicImmutableFunctionView&) MGIS_NOEXCEPT;

  /*!
   * \brief compute the Cauchy stress from the first Piola-Kirchhoff stress and
   * the deformation gradient of a function
   *
   * Both `pk1` and `F` are interpreted as unsymmetric tensors, which requires
   * 3, 5 or 9 components, corresponding to unsymmetric tensors of dimension 1,
   * 2 or 3; `pk1` and `F` must have the same number of components.
   *
   * \param[in] ctx: execution context used to report the errors
   * \param[in] pk1: view on the function whose values are the first
   * Piola-Kirchhoff stresses
   * \param[in] F: view on the function whose values are the deformation
   * gradients
   * \return the function holding the Cauchy stresses or an empty optional if
   * the number of components of `pk1` is not compatible with an unsymmetric
   * tensor, in which case an error is reported through the error handler
   */
  MGIS_EXPORT std::optional<Function<BasicLinearSpace>>
  computeCauchyStressFromFirstPiolaKirchhoffStress(
      Context&,
      const BasicImmutableFunctionView&,
      const BasicImmutableFunctionView&) MGIS_NOEXCEPT;

  /*!
   * \brief compute the hydrostatic pressure of the Cauchy stress obtained from
   * the first Piola-Kirchhoff stress and the deformation gradient of a function
   *
   * The Cauchy stress is first computed from `pk1` and `F` and its hydrostatic
   * pressure, defined as `trace(s) / 3`, is then evaluated. Both `pk1` and `F`
   * are interpreted as unsymmetric tensors, which requires 3, 5 or 9
   * components, corresponding to unsymmetric tensors of dimension 1, 2 or 3;
   * `pk1` and `F` must have the same number of components.
   *
   * \param[in] ctx: execution context used to report the errors
   * \param[in] pk1: view on the function whose values are the first
   * Piola-Kirchhoff stresses
   * \param[in] F: view on the function whose values are the deformation
   * gradients
   * \return the function holding the hydrostatic pressures or an empty optional
   * if the number of components of `pk1` is not compatible with an unsymmetric
   * tensor, in which case an error is reported through the error handler
   */
  MGIS_EXPORT std::optional<Function<BasicLinearSpace>>
  computeCauchyStressHydrostraticPressureFromFirstPiolaKirchhoffStress(
      Context&,
      const BasicImmutableFunctionView&,
      const BasicImmutableFunctionView&) MGIS_NOEXCEPT;

  /*!
   * \brief compute the von Mises stress of the Cauchy stress obtained from the
   * first Piola-Kirchhoff stress and the deformation gradient of a function
   *
   * The Cauchy stress is first computed from `pk1` and `F` and its von Mises
   * stress is then evaluated. Both `pk1` and `F` are interpreted as unsymmetric
   * tensors, which requires 3, 5 or 9 components, corresponding to unsymmetric
   * tensors of dimension 1, 2 or 3; `pk1` and `F` must have the same number of
   * components.
   *
   * \param[in] ctx: execution context used to report the errors
   * \param[in] pk1: view on the function whose values are the first
   * Piola-Kirchhoff stresses
   * \param[in] F: view on the function whose values are the deformation
   * gradients
   * \return the function holding the von Mises stresses or an empty optional if
   * the number of components of `pk1` is not compatible with an unsymmetric
   * tensor, in which case an error is reported through the error handler
   */
  MGIS_EXPORT std::optional<Function<BasicLinearSpace>>
  computeCauchyStressVonMisesStressFromFirstPiolaKirchhoffStress(
      Context&,
      const BasicImmutableFunctionView&,
      const BasicImmutableFunctionView&) MGIS_NOEXCEPT;

  /*!
   * \brief compute the eigen values of the Cauchy stress obtained from the
   * first Piola-Kirchhoff stress and the deformation gradient of a function
   *
   * The Cauchy stress is first computed from `pk1` and `F` and its eigen values
   * are then evaluated. Both `pk1` and `F` are interpreted as unsymmetric
   * tensors, which requires 3, 5 or 9 components, corresponding to unsymmetric
   * tensors of dimension 1, 2 or 3; `pk1` and `F` must have the same number of
   * components.
   *
   * \param[in] ctx: execution context used to report the errors
   * \param[in] pk1: view on the function whose values are the first
   * Piola-Kirchhoff stresses
   * \param[in] F: view on the function whose values are the deformation
   * gradients
   * \return the function holding the eigen values or an empty optional if the
   * number of components of `pk1` is not compatible with an unsymmetric tensor,
   * in which case an error is reported through the error handler
   */
  MGIS_EXPORT std::optional<Function<BasicLinearSpace>>
  computeCauchyStressEigenValuesFromFirstPiolaKirchhoffStress(
      Context&,
      const BasicImmutableFunctionView&,
      const BasicImmutableFunctionView&) MGIS_NOEXCEPT;

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TFEL_MECHANICALOPERATIONS_HXX */
