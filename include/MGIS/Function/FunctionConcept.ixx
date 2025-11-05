/*!
 * \file   MGIS/Function/FunctionConcept.ixx
 * \brief
 * \author Thomas Helfer
 * \date   14/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_FUNCTIONCONCEPT_IXX
#define LIB_MGIS_FUNCTION_FUNCTIONCONCEPT_IXX

#include "MGIS/Context.hxx"
#include "MGIS/Function/Algorithms.hxx"

namespace mgis::function::internals {

  template <FunctionConcept FunctionType>
  constexpr decltype(auto) disambiguateGetSpace(const FunctionType& f)  //
      requires(!EvaluatorConcept<FunctionType>) {
    return getSpace(f);
  }  // end of disambiguateGetSpace

  template <FunctionConcept FunctionType>
  constexpr mgis::size_type disambiguateGetNumberOfComponents(
      const FunctionType& f)  //
      requires(!EvaluatorConcept<FunctionType>) {
    return getNumberOfComponents(f);
  }  // end of disambiguateGetNumberOfComponents

}  // end of namespace mgis::function::internals

namespace mgis::function {

  template <EvaluatorConcept EvaluatorType, FunctionConcept FunctionType>
  bool operator|(EvaluatorType e, FunctionType& f) requires(
      std::same_as<evaluator_space<EvaluatorType>,
                   function_space<FunctionType>>) {
    Context ctx;
    return assign(ctx, f, e);
  }  // end of operator |

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_FUNCTIONCONCEPT_IXX */
