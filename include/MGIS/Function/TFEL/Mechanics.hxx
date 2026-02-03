/*!
 * \file   MGIS/Function/TFEL/Mechanics.hxx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef MGIS_HAVE_TFEL
#error "TFEL is required to use this header"
#endif /* MGIS_HAVE_TFEL */

#ifndef LIB_MGIS_FUNCTION_TFEL_MECHANICS_HXX
#define LIB_MGIS_FUNCTION_TFEL_MECHANICS_HXX

#include <type_traits>
#include "TFEL/Material/FiniteStrainBehaviourTangentOperator.hxx"
#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/TFEL/Tensors.hxx"

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

#include "MGIS/Function/TFEL/Mechanics.ixx"

#endif /* LIB_MGIS_FUNCTION_MECHANICS_HXX */
