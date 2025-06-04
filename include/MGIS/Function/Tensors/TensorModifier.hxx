/*!
 * \file   MGIS/Function/Tensors/TensorModifier.hxx
 * \brief
 * \author Thomas Helfer
 * \date   13/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_TENSORMODIFIER_HXX
#define LIB_MGIS_FUNCTION_TENSORMODIFIER_HXX

#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/EvaluatorModifierBase.hxx"
#include "MGIS/Function/Tensors/TensorConcept.hxx"

namespace mgis::function {

  /*
   * \brief a variable stating if the requirements of the
   * TensorView class are met.
   *
   * \tparam TensorType: tensorial object to be exposed
   * \tparam EvaluatorType: evaluator of the stress
   */
  template <TensorConcept TensorType, EvaluatorConcept EvaluatorType>
  inline constexpr bool areTensorModifierRequirementsSatisfied =
      ((isEvaluatorResultTypeMappable<EvaluatorType>)&&  //
       (number_of_components<EvaluatorType> != dynamic_extent
            ? compile_time_size<TensorType> ==
                  number_of_components<EvaluatorType>
            : true));

  /*!
   * \brief a base class for evaluators exposing tensorial object
   * \tparam TensorType: tensorial object to be exposed
   * \tparam EvaluatorType: evaluator of the stress
   */
  template <TensorConcept TensorType, EvaluatorConcept EvaluatorType>
  requires(areTensorModifierRequirementsSatisfied<TensorType,
                                                  EvaluatorType>)  //
      struct TensorModifier
      : EvaluatorModifierBase<TensorModifier<TensorType, EvaluatorType>,
                              EvaluatorType> {
    // inheriting constructor
    using EvaluatorModifierBase<TensorModifier,
                                EvaluatorType>::EvaluatorModifierBase;
    //! \return the size of the tensorial object to be exposed
    static constexpr size_type getNumberOfComponents() noexcept;
    /*!
     * \brief apply the modifier
     * \param[in] values: values to be modified
     */
    constexpr auto apply(const evaluator_result<EvaluatorType>&) const;
  };

}  // end of namespace mgis::function

#include "MGIS/Function/Tensors/TensorModifier.ixx"

#endif /* LIB_MGIS_FUNCTION_TENSORMODIFIER_HXX */
