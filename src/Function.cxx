/*!
 * \file   src/Function.cxx
 * \brief
 * \author Thomas Helfer
 * \date   7/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/Function.hxx"

namespace mgis::function {

  static_assert(EvaluatorConcept<FunctionEvaluator<BasicLinearSpace>>);
  static_assert(!FunctionConcept<FunctionEvaluator<BasicLinearSpace>>);

  static_assert(EvaluatorConcept<FunctionView<BasicLinearSpace>>);
  static_assert(FunctionConcept<FunctionView<BasicLinearSpace>>);
  static_assert(std::same_as<function_result<FunctionView<BasicLinearSpace>>,
                             std::span<real>>);
  static_assert(
      std::same_as<function_result<FunctionView<BasicLinearSpace,
                                                FunctionDataLayoutDescription{
                                                    .data_size = 2}>>,
                   std::span<real, 2>>);
  static_assert(
      std::same_as<function_result<FunctionView<BasicLinearSpace,
                                                FunctionDataLayoutDescription{
                                                    .data_size = 1}>>,
                   real&>);

  static_assert(FunctionConcept<Function<BasicLinearSpace>>);
  // This shall not work as Function is not a lightweight obect
  static_assert(!EvaluatorConcept<Function<BasicLinearSpace>>);

  static_assert(EvaluatorConcept<FixedSizeView<Function<BasicLinearSpace>, 9>>);
  static_assert(FunctionConcept<FixedSizeView<Function<BasicLinearSpace>, 9>>);

}  // end of namespace mgis::function
