/*!
 * \file   UniformEvaluator.cxx
 * \brief
 * \author Thomas Helfer
 * \date   06/10/2025
 */

#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/UniformEvaluator.hxx"

namespace mgis::function {

  static_assert(EvaluatorConcept<UniformEvaluator<BasicLinearSpace, 1>>);
  static_assert(EvaluatorConcept<UniformEvaluator<BasicLinearSpace, 6>>);

}  // end of namespace mgis::function
