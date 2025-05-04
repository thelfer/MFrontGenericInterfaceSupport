/*!
 * \file   Evaluators.cxx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/Evaluators.hxx"

namespace mgis::function {

  static_assert(EvaluatorConcept<FixedSizedEvaluator<BasicLinearSpace, 9>>);

}  // end of namespace mgis::function