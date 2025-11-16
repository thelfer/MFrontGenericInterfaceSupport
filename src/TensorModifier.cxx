/*!
 * \file   TensorModifier.cxx
 * \brief
 * \author Thomas Helfer
 * \date   01/09/2025
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
#include "MGIS/Function/TFEL/Tensors.hxx"
#include "MGIS/Function/TFEL/TensorModifier.hxx"

namespace mgis::function {

  static_assert(
      EvaluatorConcept<TensorModifier<tfel::math::stensor<3u, mgis::real>,
                                      FunctionView<BasicLinearSpace>>>);
  static_assert(
      !FunctionConcept<TensorModifier<tfel::math::stensor<3u, mgis::real>,
                                      FunctionView<BasicLinearSpace>>>);

  static_assert(
      number_of_components<TensorModifier<tfel::math::stensor<3u, mgis::real>,
                                          FunctionView<BasicLinearSpace>>> ==
      6);

}  // end of namespace mgis::function
