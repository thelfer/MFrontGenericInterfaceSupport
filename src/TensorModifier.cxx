/*!
 * \file   TensorModifier.cxx
 * \brief
 * \author Thomas Helfer
 * \date   01/09/2025
 */

#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/Tensors.hxx"
#include "MGIS/Function/Tensors/TensorModifier.hxx"

namespace mgis::function {

#ifdef MGIS_HAVE_TFEL

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

#endif /* MGIS_HAVE_TFEL */

}  // end of namespace mgis::function
