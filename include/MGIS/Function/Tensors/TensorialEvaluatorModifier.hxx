/*!
 * \file   MGIS/Function/TensorialEvaluatorModifier.hxx
 * \brief  This file declares the TensorialEvaluatorModifier class and the
 * transform function, as well as the as_fsarray, as_tvector,
 * as_tmatrix, as_stensor and as_tensor modifiers.
 * \author Thomas Helfer
 * \date   09/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_TENSORIALEVALUATORMODIFIER_HXX
#define LIB_MGIS_FUNCTION_TENSORIALEVALUATORMODIFIER_HXX

#include <span>
#include "MGIS/Function/EvaluatorModifierBase.hxx"
#include "MGIS/Function/Tensors/TensorialObject.hxx"

namespace mgis::function::internals {

  /*
   * \brief a variable stating if the requirements of the
   * TensorialEvaluatorModifier class are met.
   *
   * \tparam TensorType: tensorial object to be exposed
   * \tparam EvaluatorType: evaluator of the stress
   */
  template <TensorialObjectConcept TensorType, EvaluatorConcept EvaluatorType>
  inline constexpr bool areTensorialEvaluatorModifierRequirementsSatisfied =
      ((std::same_as<evaluator_result<EvaluatorType>, std::span<const real>>) ||
       (std::same_as<evaluator_result<EvaluatorType>,
                     std::span<const real, compile_time_size<TensorType>>>)) &&
      (number_of_components<EvaluatorType> != dynamic_extent
           ? compile_time_size<TensorType> ==
                 number_of_components<EvaluatorType>
           : true);

  /*!
   * \brief a base class for evaluators exposing tensorial object
   * \tparam TensorType: tensorial object to be exposed
   * \tparam EvaluatorType: evaluator of the stress
   */
  template <TensorialObjectConcept TensorType, EvaluatorConcept EvaluatorType>
  requires(
      areTensorialEvaluatorModifierRequirementsSatisfied<TensorType,
                                                         EvaluatorType>)  //
      struct TensorialEvaluatorModifier
      : EvaluatorModifierBase<
            TensorialEvaluatorModifier<TensorType, EvaluatorType>,
            EvaluatorType> {
    // inheriting constructor
    using EvaluatorModifierBase<TensorialEvaluatorModifier,
                                EvaluatorType>::EvaluatorModifierBase;
    //! \return the size of the tensorial object to be exposed
    static constexpr size_type getNumberOfComponents() noexcept;
    /*!
     * \brief apply the modifier
     * \param[in] values: values to be modified
     */
    auto apply(const evaluator_result<EvaluatorType>&) const;
  };

  template <TensorialObjectConcept TensorType>
  struct tensor_modifier {
    template <typename EvaluatorType>
    auto operator()(EvaluatorType&&) const
        requires((EvaluatorConcept<std::decay_t<EvaluatorType>>)&&(
            areTensorialEvaluatorModifierRequirementsSatisfied<TensorType,
                                                               EvaluatorType>));
  };

}  // end of namespace mgis::function::internals

namespace mgis::function {

  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      inline constexpr auto as_stensor =
          internals::tensor_modifier<tfel::math::stensor<N, real>>{};

  template <unsigned short N>
  requires((N == 1) || (N == 2) || (N == 3))  //
      inline constexpr auto as_tensor =
          internals::tensor_modifier<tfel::math::tensor<N, real>>{};

  template <unsigned short N>
  inline constexpr auto as_fsarray =
      internals::tensor_modifier<tfel::math::fsarray<N, real>>{};

  template <unsigned short N>
  inline constexpr auto as_tvector =
      internals::tensor_modifier<tfel::math::tvector<N, real>>{};

  template <unsigned short N, unsigned short M>
  inline constexpr auto as_tmatrix =
      internals::tensor_modifier<tfel::math::tmatrix<N, M, real>>{};

}  // end of namespace mgis::function

#include "MGIS/Function/Tensors/TensorialEvaluatorModifier.ixx"

#endif /* LIB_MGIS_FUNCTION_TENSORIALEVALUATORMODIFIER_HXX */
