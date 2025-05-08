/*!
 * \file   src/Function.cxx
 * \brief
 * \author Thomas Helfer
 * \date   7/05/2025
 */

#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/Function.hxx"

namespace mgis::function {

  static_assert(EvaluatorConcept<ImmutableFunctionView<BasicLinearSpace>>);

#pragma message("This shall not work as Function is not a lightweight obect")
  static_assert(EvaluatorConcept<Function<BasicLinearSpace>>);

}  // end of namespace mgis::function
