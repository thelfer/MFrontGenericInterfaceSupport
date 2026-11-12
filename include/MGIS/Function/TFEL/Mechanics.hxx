/*!
 * \file   MGIS/Function/TFEL/Mechanics.hxx
 * \brief  File which defines the mechanical operations which can be applied to
 * evaluators returning tensorial values (`hydrostatic_stress`, `vmis`,
 * `principal_stress`, `from_pk1_to_cauchy`, ...) and the
 * `convert_finite_strain_stiffness` function converting finite strain tangent
 * operators from one kind to another.
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
#include "MGIS/Function/TFEL/General.hxx"
#include "MGIS/Function/TFEL/Tensors.hxx"

namespace mgis::function::internals {

  /*!
   * \brief callable object computing the hydrostatic stress, i.e. the mean
   * normal stress, of a symmetric tensor, defined as `trace(s) / 3`
   */
  struct HydrostaticStressOperator {
    /*!
     * \brief compute the hydrostatic stress of a symmetric tensor
     * \param[in] s: symmetric tensor considered
     * \return `trace(s) / 3`
     */
    MGIS_HOST_DEVICE [[nodiscard]] constexpr auto operator()(
        const tfel::math::StensorConcept auto& s) const MGIS_NOEXCEPT {
      return tfel::math::trace(s) / 3;
    }
  };

  /*!
   * \brief callable object computing the von Mises stress of a symmetric tensor
   * interpreted as a stress
   */
  struct VonMisesStressOperator {
    /*!
     * \brief compute the von Mises stress of a symmetric tensor
     * \param[in] s: symmetric tensor interpreted as a stress
     * \return the von Mises stress, computed by `tfel::math::sigmaeq`
     */
    MGIS_HOST_DEVICE [[nodiscard]] constexpr auto operator()(
        const tfel::math::StensorConcept auto& s) const MGIS_NOEXCEPT {
      return tfel::math::sigmaeq(s);
    }
  };

  /*!
   * \brief callable object converting a first Piola-Kirchhoff stress into a
   * Cauchy stress given the deformation gradient
   */
  struct FirstPiolaKirchhoffStressToCauchyStressOperator {
    /*!
     * \brief convert a first Piola-Kirchhoff stress into a Cauchy stress
     * \tparam FirstPiolaKirchhoffStressType: type of the first Piola-Kirchhoff
     * stress, an unsymmetric tensor
     * \tparam DeformationGradientType: type of the deformation gradient, an
     * unsymmetric tensor
     * \param[in] pk1: first Piola-Kirchhoff stress
     * \param[in] F: deformation gradient
     * \return the corresponding Cauchy stress, a symmetric tensor
     * \note only available if the two tensorial arguments have the same space
     * dimension
     */
    template <tfel::math::TensorConcept FirstPiolaKirchhoffStressType,
              tfel::math::TensorConcept DeformationGradientType>
    MGIS_HOST_DEVICE [[nodiscard]] constexpr auto operator()(
        const FirstPiolaKirchhoffStressType& pk1,
        const DeformationGradientType& F) const MGIS_NOEXCEPT
      requires(tfel::math::getSpaceDimension<FirstPiolaKirchhoffStressType>() ==
               tfel::math::getSpaceDimension<DeformationGradientType>())
    {
      return tfel::math::convertFirstPiolaKirchhoffStressToCauchyStress(pk1, F);
    }
  };

}  // end of namespace mgis::function::internals

namespace mgis::function {

  /*!
   * \brief modifier computing the hydrostatic stress (`trace(s) / 3`) of the
   * symmetric tensorial values computed by an evaluator
   */
  inline constexpr auto hydrostatic_stress =
      internals::unary_operation_modifier2(
          internals::HydrostaticStressOperator{});

  /*!
   * \brief modifier computing the von Mises stress of the symmetric tensorial
   * values computed by an evaluator
   */
  inline constexpr auto vmis =
      internals::unary_operation_modifier2(internals::VonMisesStressOperator{});

  //! \brief alias for `vmis`, modifier computing the von Mises stress
  inline constexpr auto von_mises_stress = vmis;

  /*!
   * \brief modifier computing the principal stresses, i.e. the eigen values, of
   * the symmetric tensorial values computed by an evaluator
   *
   * \tparam esolver: eigen solver used to compute the eigen values, which
   * defaults to `tfel::math::stensor_common::TFELEIGENSOLVER`
   */
  template <tfel::math::stensor_common::EigenSolver esolver =
                tfel::math::stensor_common::TFELEIGENSOLVER>
  inline constexpr auto principal_stress = eigen_values<esolver>;

  /*!
   * \brief modifier converting the first Piola-Kirchhoff stresses computed by
   * an evaluator into Cauchy stresses, given an evaluator of the deformation
   * gradient
   */
  inline constexpr auto from_pk1_to_cauchy =
      internals::binary_operation_modifier2(
          internals::FirstPiolaKirchhoffStressToCauchyStressOperator{});

  //! \brief conversion of finite strain tangent operators
  using FiniteStrainStiffnessKind =
      tfel::material::FiniteStrainBehaviourTangentOperatorBase::Flag;

  /*!
   * \brief create a modifier converting a finite strain tangent operator
   * computed by an evaluator into another kind of tangent operator
   *
   * The conversion is performed by
   * `tfel::material::convert<ResultFlag, SourceFlag, N>` and generally requires
   * the two deformation gradients and the Cauchy stress at the integration
   * point considered. The returned object is a modifier to be applied, using
   * `operator|`, to an evaluator returning a tangent operator of kind
   * `SourceFlag`.
   *
   * \tparam ResultFlag: kind of the tangent operator to be returned
   * \tparam SourceFlag: kind of the tangent operator taken as input
   * \tparam DeformationGradientEvaluatorType0: type of the evaluator providing
   * the first deformation gradient
   * \tparam DeformationGradientEvaluatorType1: type of the evaluator providing
   * the second deformation gradient
   * \tparam CauchyStressEvaluatorType: type of the evaluator providing the
   * Cauchy stress
   * \param[in] F0: evaluator used to compute the first deformation gradient
   * \param[in] F1: evaluator used to compute the second deformation gradient
   * \param[in] s: evaluator used to compute the Cauchy stress
   * \return a modifier converting the tangent operators it is applied to
   */
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
