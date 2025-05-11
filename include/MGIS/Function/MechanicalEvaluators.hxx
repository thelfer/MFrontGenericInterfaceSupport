/*!
 * \file   MGIS/Function/MechanicalEvaluators.hxx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#ifndef LIB_FUNCTION_MECHANICALEVALUATORS_HXX
#define LIB_FUNCTION_MECHANICALEVALUATORS_HXX

#include "MGIS/Function/Evaluator.hxx"

#ifdef MGIS_HAVE_TFEL

#include "MGIS/Function/TransformEvaluatorModifier.hxx"
#include "MGIS/Function/Tensors/TensorialObject.hxx"
#include "MGIS/Function/Tensors/SymmetricTensorModifiers.hxx"

namespace mgis::function{

  inline constexpr auto hydrostatic_stress =
      internals::transform_modifier2<decltype(
          [](const tfel::math::StensorConcept auto& s) {
            return tfel::math::trace(s) / 3;
          })>{};

  inline constexpr auto vmis = internals::transform_modifier2<decltype(
      [](const tfel::math::StensorConcept auto& s) {
        return tfel::math::sigmaeq(s);
      })>{};

  template <tfel::math::stensor_common::EigenSolver esolver =
                tfel::math::stensor_common::TFELEIGENSOLVER>
  inline constexpr auto principal_stress = eigen_values<esolver>;

}

#include "MGIS/Function/MechanicalEvaluators/CauchyStressFromFirstPiolaKirchhoffStressEvaluator.hxx"

#endif MGIS_HAVE_TFEL

#endif /* LIB_FUNCTION_MECHANICALEVALUATORS_HXX */
