/*!
 * \file   MGIS/Function/Evaluator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   29/04/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_EVALUATOR_HXX
#define LIB_MGIS_FUNCTION_EVALUATOR_HXX

#include <type_traits>
#include "MGIS/Config.hxx"
#include "MGIS/AbstractErrorHandler.hxx"
#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/CompileTimeSize.hxx"
#include "MGIS/Function/EvaluatorConcept.hxx"
#include "MGIS/Function/EvaluatorModifierConcept.hxx"

namespace mgis::function {

  /*!
   * \brief check if the given evaluators shares the same space
   *
   * \param[in] ctx: context
   * \param[in] e1: first evaluator
   * \param[in] e2: second evaluator
   */
  [[nodiscard]] constexpr bool checkMatchingSpaces(
      AbstractErrorHandler&,
      const EvaluatorConcept auto&,
      const EvaluatorConcept auto&);

}  // end of namespace mgis::function

#include "MGIS/Function/Evaluator.ixx"
#include "MGIS/Function/UnaryOperation.hxx"
#include "MGIS/Function/BinaryOperation.hxx"

#endif /* LIB_MGIS_FUNCTION_EVALUATOR_HXX */
