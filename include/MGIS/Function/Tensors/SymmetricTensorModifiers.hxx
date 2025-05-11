/*!
 * \file   MGIS/Function/Tensors/SymmetricTensorModifiers.hxx
 * \brief  This headers declare
 * \author Thomas Helfer
 * \date   10/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_TENSORS_SYMMETRICTENSORMODIFIERS_HXX
#define LIB_MGIS_FUNCTION_TENSORS_SYMMETRICTENSORMODIFIERS_HXX

#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/TransformEvaluatorModifier.hxx"
#include "MGIS/Function/Tensors/TensorialObject.hxx"

namespace mgis::function {

  inline constexpr auto trace = internals::transform_modifier2<decltype(
      [](const tfel::math::StensorConcept auto& s) {
        return tfel::math::trace(s);
      })>{};

  inline constexpr auto det = internals::transform_modifier2<decltype(
      [](const tfel::math::StensorConcept auto& s) {
        return tfel::math::det(s);
      })>{};

  template <tfel::math::stensor_common::EigenSolver esolver =
                tfel::math::stensor_common::TFELEIGENSOLVER>
  inline constexpr auto eigen_values = internals::transform_modifier2<decltype(
      [](const tfel::math::StensorConcept auto& s) {
        auto s2 = eval(s);
        return s2.template computeEigenValues<esolver>();
      })>{};

} // end of mgis::function

#endif /* LIB_MGIS_FUNCTION_TENSORS_SYMMETRICTENSORMODIFIERS_HXX */
