/*!
 * \file   MGIS/Function/EvaluatorModifierConcept.ixx
 * \brief
 * \author Thomas Helfer
 * \date   05/06/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_EVALUATORMODIFIERCONCEPT_IXX
#define LIB_MGIS_FUNCTION_EVALUATORMODIFIERCONCEPT_IXX

namespace mgis::function {

  template <EvaluatorConcept EvaluatorType,
            EvaluatorModifierConcept EvaluatorModifierType>
  constexpr auto operator|(const EvaluatorType& e, EvaluatorModifierType m) {
    return m(e);
  }  // end of operator |

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_EVALUATORMODIFIERCONCEPT_IXX */
