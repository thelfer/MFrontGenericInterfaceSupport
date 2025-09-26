/*!
 * \file   MGIS/Function/Mechanics.hxx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_MECHANICS_HXX
#define LIB_MGIS_FUNCTION_MECHANICS_HXX

#include <type_traits>
#ifdef MGIS_HAVE_TFEL
#include "TFEL/Material/FiniteStrainBehaviourTangentOperator.hxx"
#endif /* MGIS_HAVE_TFEL */

#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/Tensors.hxx"

#ifdef MGIS_HAVE_TFEL

namespace mgis::function {

  inline constexpr auto hydrostatic_stress =
      internals::unary_operation_modifier2(
          [](const tfel::math::StensorConcept auto& s) {
            return tfel::math::trace(s) / 3;
          });

  inline constexpr auto vmis = internals::unary_operation_modifier2(
      [](const tfel::math::StensorConcept auto& s) {
        return tfel::math::sigmaeq(s);
      });

  inline constexpr auto von_mises_stress = vmis;

  template <tfel::math::stensor_common::EigenSolver esolver =
                tfel::math::stensor_common::TFELEIGENSOLVER>
  inline constexpr auto principal_stress = eigen_values<esolver>;

  inline constexpr auto from_pk1_to_cauchy =
      internals::binary_operation_modifier2(
          []<tfel::math::TensorConcept TensorType>(const TensorType& pk1,
                                                   const TensorType& F) {
            return tfel::math::convertFirstPiolaKirchhoffStressToCauchyStress(
                pk1, F);
          });

  // convertion of finite strain tangent operators
  using FiniteStrainStiffnessKind =
      tfel::material::FiniteStrainBehaviourTangentOperatorBase::Flag;

  template <FiniteStrainStiffnessKind ResultFlag,
            FiniteStrainStiffnessKind SourceFlag,
            TensorEvaluatorConcept DeformationGradientEvaluatorType0,
            TensorEvaluatorConcept DeformationGradientEvaluatorType1,
            StensorEvaluatorConcept CauchyStressEvaluatorType>
  constexpr auto convert_finite_strain_stiffness(
      const DeformationGradientEvaluatorType0&,
      const DeformationGradientEvaluatorType1&,
      const CauchyStressEvaluatorType&);

}  // end of namespace mgis::function

#endif /* MGIS_HAVE_TFEL */

#include "MGIS/Function/Mechanics.ixx"

#endif /* LIB_MGIS_FUNCTION_MECHANICS_HXX */
