/*!
 * \file   Evaluators.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#include "MGIS/QuadratureFunction/Evaluators.hxx"

namespace mgis::quadrature_function {

  static_assert(EvaluatorConcept<FixedSizedEvaluator<9>>);

}  // end of namespace mgis::quadrature_function