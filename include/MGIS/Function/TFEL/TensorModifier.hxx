/*!
 * \file   MGIS/Function/TFEL/TensorModifier.hxx
 * \brief
 * \author Thomas Helfer
 * \date   13/05/2025
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

#ifndef LIB_MGIS_FUNCTION_TFEL_TENSORMODIFIER_HXX
#define LIB_MGIS_FUNCTION_TFEL_TENSORMODIFIER_HXX

#include "MGIS/Function/Evaluator.hxx"
#include "MGIS/Function/EvaluatorModifierBase.hxx"
#include "MGIS/Function/TFEL/TensorConcept.hxx"

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
    constexpr size_type getNumberOfComponents() noexcept;
    /*!
     * \brief apply the modifier
     * \param[in] values: values to be modified
     */
    constexpr auto apply(const evaluator_result<EvaluatorType>&) const;
  };

  template <TensorConcept TensorType, EvaluatorConcept EvaluatorType>
  constexpr mgis::size_type getNumberOfComponents(
      const TensorModifier<TensorType, EvaluatorType>&) noexcept;

}  // end of namespace mgis::function

#include "MGIS/Function/TFEL/TensorModifier.ixx"

#endif /* LIB_MGIS_FUNCTION_TFEL_TENSORMODIFIER_HXX */
