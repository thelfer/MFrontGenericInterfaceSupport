/*!
 * \file   src/Function.cxx
 * \brief
 * \author Thomas Helfer
 * \date   7/05/2025
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
      std::same_as<
          function_result<FunctionView<BasicLinearSpace, {.data_size = 2}>>,
          std::span<real, 2>>);
  static_assert(
      std::same_as<
          function_result<FunctionView<BasicLinearSpace, {.data_size = 1}>>,
          real&>);

  // This shall not work as Function is not a lightweight obect
  static_assert(!EvaluatorConcept<Function<BasicLinearSpace>>);

  // This shall not work as Function is not a lightweight obect
  static_assert(FunctionConcept<Function<BasicLinearSpace>>);

}  // end of namespace mgis::function
