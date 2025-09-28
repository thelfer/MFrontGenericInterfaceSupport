/*!
 * \file   TensorModifier.cxx
 * \brief
 * \author Thomas Helfer
 * \date   01/09/2025
 */

#ifdef MGIS_HAVE_TFEL

#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/Tensors.hxx"
#include "MGIS/Function/Tensors/TensorModifier.hxx"

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

#endif /* MGIS_HAVE_TFEL */
