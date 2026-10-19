/*!
 * \file   MGIS/Function/TFEL/TensorModifier.ixx
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

#ifndef LIB_MGIS_FUNCTION_TFEL_TENSORMODIFIER_IXX
#define LIB_MGIS_FUNCTION_TFEL_TENSORMODIFIER_IXX

namespace mgis::function {

  template <TensorConcept TensorType, EvaluatorConcept EvaluatorType>
  requires(areTensorModifierRequirementsSatisfied<TensorType,
                                                  EvaluatorType>)  //
      constexpr size_type
      TensorModifier<TensorType,
                     EvaluatorType>::getNumberOfComponents() noexcept {
    return compile_time_size<TensorType>;
  }

  //! \brief perform consistency checks
  template <TensorConcept TensorType, EvaluatorConcept EvaluatorType>
  requires(areTensorModifierRequirementsSatisfied<TensorType,
                                                  EvaluatorType>)  //
      constexpr bool TensorModifier<TensorType, EvaluatorType>::check(
          AbstractErrorHandler& eh) const {
    const auto c =
        EvaluatorModifierBase<TensorModifier<TensorType, EvaluatorType>,
                              EvaluatorType>::check(eh);
    if (!c) {
      return false;
    }
    const auto nc =
        internals::disambiguateGetNumberOfComponents(this->evaluator);
    if (nc != compile_time_size<TensorType>) {
      return eh.registerErrorMessage("invalid number of components");
    }
    return true;
  }

  template <TensorConcept TensorType, EvaluatorConcept EvaluatorType>
  requires(areTensorModifierRequirementsSatisfied<TensorType,
                                                  EvaluatorType>)  //
      constexpr auto TensorModifier<TensorType, EvaluatorType>::apply(
          const evaluator_result<EvaluatorType>& values) const {
    if constexpr (DataViewConcept<evaluator_result<EvaluatorType>>) {
      return tfel::math::map<const TensorType>(values.data());
    } else {
      return TensorType(tfel::math::map<const TensorType>(values.data()));
    }
  }  // end of apply

  template <TensorConcept TensorType, EvaluatorConcept EvaluatorType>
  constexpr mgis::size_type getNumberOfComponents(
      const TensorModifier<TensorType, EvaluatorType>& e) noexcept {
    return e.getNumberOfComponents();
  }

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TFEL_TENSORMODIFIER_IXX */
