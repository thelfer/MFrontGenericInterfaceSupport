/*!
 * \file   src/Quantity.cxx
 * \brief
 * \author Thomas Helfer
 * \date   16/08/2026
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
#include "MGIS/Function/TFEL/Quantity.hxx"

namespace mgis::function {

  static_assert(EvaluatorConcept<QuantityView<Function<BasicLinearSpace>,
                                              tfel::math::unit::Time>>);
  static_assert(FunctionConcept<QuantityView<Function<BasicLinearSpace>,
                                             tfel::math::unit::Time>>);
  static_assert(number_of_components<QuantityView<Function<BasicLinearSpace>,
                                                  tfel::math::unit::Stress>> ==
                1);

  static_assert(
      EvaluatorConcept<QuantityModifier<FunctionView<BasicLinearSpace>,
                                        tfel::math::unit::Time>>);
  static_assert(
      !FunctionConcept<QuantityModifier<FunctionView<BasicLinearSpace>,
                                        tfel::math::unit::Time>>);

  static_assert(
      number_of_components<QuantityModifier<FunctionView<BasicLinearSpace>,
                                            tfel::math::unit::Stress>> == 1);

}  // end of namespace mgis::function
