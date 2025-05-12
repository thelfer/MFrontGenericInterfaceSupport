/*!
 * \file   MGIS/Function/Tensors.hxx
 * \brief
 * \author Thomas Helfer
 * \date   11/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_TENSORS_HXX
#define LIB_MGIS_FUNCTION_TENSORS_HXX

#include "MGIS/Function/Evaluator.hxx"

#ifdef MGIS_HAVE_TFEL

#include "MGIS/Function/Tensors/TensorialObject.hxx"

namespace mgis::function {

  inline constexpr auto trace = internals::unary_operation_modifier2(
      []<typename TensorType>(const TensorType& t) requires(
          (tfel::math::TensorConcept<TensorType>) ||
          (tfel::math::StensorConcept<TensorType>)) {
        return tfel::math::trace(t);
      });

  inline constexpr auto det = internals::unary_operation_modifier2(
      []<typename TensorType>(const TensorType& t) requires(
          (tfel::math::TensorConcept<TensorType>) ||
          (tfel::math::StensorConcept<TensorType>)) {
        return tfel::math::det(t);
      });

  template <tfel::math::stensor_common::EigenSolver esolver =
                tfel::math::stensor_common::TFELEIGENSOLVER>
  inline constexpr auto eigen_values = internals::unary_operation_modifier2(
      [](const tfel::math::StensorConcept auto& s) {
        auto s2 = eval(s);
        return s2.template computeEigenValues<esolver>();
      });

  inline constexpr auto rotate = internals::binary_operation_modifier2(
      []<typename TensorType>(const TensorType& pk1,
                              const tfel::math::tmatrix<3, 3>& R)  //
      requires((tfel::math::TensorConcept<TensorType>) ||
               (tfel::math::StensorConcept<TensorType>)) {
        return tfel::math::change_basis(pk1, R);
      });

  inline constexpr auto rotate_backwards =
      internals::binary_operation_modifier2(
          []<typename TensorType>(const TensorType& pk1,
                                  const tfel::math::tmatrix<3, 3>& R)  //
          requires((tfel::math::TensorConcept<TensorType>) ||
                   (tfel::math::StensorConcept<TensorType>)) {
            return tfel::math::change_basis(pk1, transpose(R));
          });

  template <tfel::math::stensor_common::EigenSolver esolver =
                tfel::math::stensor_common::TFELEIGENSOLVER>
  inline constexpr auto logarithmic_strain =
      internals::unary_operation_modifier2(
          [](const tfel::math::TensorConcept auto& F) {
            using namespace tfel::math;
            constexpr auto N = getSpaceDimension<decltype(F)>();
            const auto e = computeGreenLagrangeTensor(F);
            const auto [vp, m] = e.template computeEigenVectors<esolver>();
            const auto log_vp =
                map([](const real x) { return std::log1p(2 * x) / 2; }, vp);
            return stensor<N, real>::computeIsotropicFunction(e, m);
          });

}  // end of namespace mgis::function

#endif MGIS_HAVE_TFEL

#endif /* LIB_MGIS_FUNCTION_TENSORS_HXX */
