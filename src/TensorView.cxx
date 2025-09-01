/*!
 * \file   TensorView.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   01/09/2025
 */

#include "MGIS/Function/BasicLinearSpace.hxx"
#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/Tensors.hxx"
#include "MGIS/Function/Tensors/TensorView.hxx"

namespace mgis::function{

#ifdef MGIS_HAVE_TFEL

  static_assert(
      EvaluatorConcept<TensorView<Function<BasicLinearSpace>,
                                  tfel::math::stensor<3u, mgis::real>>>);
  static_assert(
      FunctionConcept<TensorView<Function<BasicLinearSpace>,
                                 tfel::math::stensor<3u, mgis::real>>>);

  static_assert(
      number_of_components<TensorView<Function<BasicLinearSpace>,
                                      tfel::math::stensor<3u, mgis::real>>> ==
      6);

#endif /* MGIS_HAVE_TFEL */

} // end of namespace mgis::function
