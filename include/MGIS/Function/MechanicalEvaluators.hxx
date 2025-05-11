/*!
 * \file   MGIS/Function/MechanicalEvaluators.hxx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_MECHANICALEVALUATORS_HXX
#define LIB_MGIS_FUNCTION_MECHANICALEVALUATORS_HXX

#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/TransformEvaluatorModifier.hxx"
#include "MGIS/Function/BinaryOperationEvaluator.hxx"
#include "MGIS/Function/Tensors.hxx"

#ifdef MGIS_HAVE_TFEL

#include "MGIS/Function/Tensors/TensorialObject.hxx"

namespace mgis::function{

  inline constexpr auto hydrostatic_stress =
      internals::transform_modifier2(
          [](const tfel::math::StensorConcept auto& s) {
            return tfel::math::trace(s) / 3;
          });

  inline constexpr auto vmis = internals::transform_modifier2(
      [](const tfel::math::StensorConcept auto& s) {
        return tfel::math::sigmaeq(s);
      });

  template <tfel::math::stensor_common::EigenSolver esolver =
                tfel::math::stensor_common::TFELEIGENSOLVER>
  inline constexpr auto principal_stress = eigen_values<esolver>;

  inline constexpr auto from_pk1_to_cauchy =
      internals::binary_operation_modifier2<decltype(
          []<tfel::math::TensorConcept TensorType>(const TensorType& pk1,
                                                   const TensorType& F) {
            return tfel::math::convertFirstPiolaKirchhoffStressToCauchyStress(
                pk1, F);
          })>{};

} // end of namespace mgis::function

#endif MGIS_HAVE_TFEL

#endif /* LIB_MGIS_FUNCTION_MECHANICALEVALUATORS_HXX */
