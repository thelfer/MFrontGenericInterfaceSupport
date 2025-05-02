/*!
 * \file   Evaluators.cxx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#include "MGIS/Function/Evaluators.hxx"

namespace mgis::function {

  static_assert(EvaluatorConcept<FixedSizedEvaluator<9>>);

}  // end of namespace mgis::function