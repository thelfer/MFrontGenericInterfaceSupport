/*!
 * \file   MGIS/Function/MechanicalOperations.hxx
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

#ifndef LIB_MGIS_FUNCTION_MECHANICALOPERATIONS_HXX
#define LIB_MGIS_FUNCTION_MECHANICALOPERATIONS_HXX

#include <optional>
#include "MGIS/Config.hxx"
#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/Function.hxx"

namespace mgis::function {

  using BasicImmutableFunctionView =
      FunctionView<BasicLinearSpace, FunctionDataLayoutDescription{}, false>;

  MGIS_EXPORT std::optional<Function<BasicLinearSpace>> computeTrace(
      Context&, const BasicImmutableFunctionView&) noexcept;

  MGIS_EXPORT std::optional<Function<BasicLinearSpace>> computeVonMisesStress(
      Context&, const BasicImmutableFunctionView&) noexcept;

  MGIS_EXPORT std::optional<Function<BasicLinearSpace>>
  computeHydrostraticPressure(Context&,
                              const BasicImmutableFunctionView&) noexcept;

  MGIS_EXPORT std::optional<Function<BasicLinearSpace>> computeEigenValues(
      Context&, const BasicImmutableFunctionView&) noexcept;

  MGIS_EXPORT std::optional<Function<BasicLinearSpace>>
  computeCauchyStressFromFirstPiolaKirchhoffStress(
      Context&,
      const BasicImmutableFunctionView&,
      const BasicImmutableFunctionView&) noexcept;

  MGIS_EXPORT std::optional<Function<BasicLinearSpace>>
  computeCauchyStressHydrostraticPressureFromFirstPiolaKirchhoffStress(
      Context&,
      const BasicImmutableFunctionView&,
      const BasicImmutableFunctionView&) noexcept;

  MGIS_EXPORT std::optional<Function<BasicLinearSpace>>
  computeCauchyStressVonMisesStressFromFirstPiolaKirchhoffStress(
      Context&,
      const BasicImmutableFunctionView&,
      const BasicImmutableFunctionView&) noexcept;

  MGIS_EXPORT std::optional<Function<BasicLinearSpace>>
  computeCauchyStressEigenValuesFromFirstPiolaKirchhoffStress(
      Context&,
      const BasicImmutableFunctionView&,
      const BasicImmutableFunctionView&) noexcept;

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_MECHANICALOPERATIONS_HXX */
