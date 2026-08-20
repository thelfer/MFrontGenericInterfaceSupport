/*!
 * \file   MGIS/Function/TFEL/QuantityModifier.ixx
 * \brief
 * \author Thomas Helfer
 * \date   16/08/2026
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_TFEL_QUANTITYMODIFIER_IXX
#define LIB_MGIS_FUNCTION_TFEL_QUANTITYMODIFIER_IXX

#include "MGIS/Raise.hxx"

namespace mgis::function {

  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr bool QuantityModifier<EvaluatorType, UnitType>::checkPreconditions(
      AbstractErrorHandler& eh, const EvaluatorType& values) {
    if (internals::disambiguateGetNumberOfComponents(values) != 1) {
      return eh.registerErrorMessage("invalid number of components");
    }
    return true;
  }  // end of checkPreconditions

  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr QuantityModifier<EvaluatorType, UnitType>::QuantityModifier(
      const EvaluatorType& e)
      : QuantityModifier(preconditions_check, e) {}  // end of QuantityModifier

  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  template <bool doPreconditionsCheck>
  constexpr QuantityModifier<EvaluatorType, UnitType>::QuantityModifier(
      const PreconditionsCheck<doPreconditionsCheck>& pcheck,
      const EvaluatorType& e)
      : PreconditionsChecker<QuantityModifier>(pcheck, e),
        evaluator(e) {}  // end of QuantityModifier

  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr bool QuantityModifier<EvaluatorType, UnitType>::check(
      AbstractErrorHandler& ctx) const {
    return checkPreconditions(ctx, this->evaluator);
  }

  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  decltype(auto) QuantityModifier<EvaluatorType, UnitType>::getSpace() const {
    return internals::disambiguateGetSpace(this->evaluator);
  }

  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr size_type
  QuantityModifier<EvaluatorType, UnitType>::getNumberOfComponents() const {
    return 1;
  }

  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr auto QuantityModifier<EvaluatorType, UnitType>::operator()(
      const element_index<Space>& i) const
      requires((internals::EvaluatorResultQuery<EvaluatorType>::b1) &&
               (isEvaluatorResultTypeMappable<EvaluatorType>)) {
    constexpr auto has_data_method =
        requires(const EvaluatorType& rf, const element_index<Space>& ri) {
      { rf.data(unsafe, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::const_qt_ref<UnitType, real>(
          *(this->evaluator.data(unsafe, i)));
    } else {
      return tfel::math::const_qt_ref<UnitType, real>(
          *(this->evaluator(i).data()));
    }
  }

  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr auto QuantityModifier<EvaluatorType, UnitType>::operator()(
      const element_workspace<Space>& wk, const element_index<Space>& i) const
      requires((internals::EvaluatorResultQuery<EvaluatorType>::b2) &&
               (isEvaluatorResultTypeMappable<EvaluatorType>)) {
    constexpr auto has_data_method =
        requires(const EvaluatorType& rf, const element_workspace<Space>& rwk,
                 const element_index<Space>& ri) {
      { rf.data(unsafe, rwk, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::const_qt_ref<UnitType, real>(
          *(this->evaluator.data(unsafe, wk, i)));
    } else {
      return tfel::math::const_qt_ref<UnitType, real>(
          *(this->evaluator(wk, i).data()));
    }
  }

  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr auto QuantityModifier<EvaluatorType, UnitType>::operator()(
      const cell_index<Space>& e, const quadrature_point_index<Space>& i) const
      requires((internals::EvaluatorResultQuery<EvaluatorType>::b3) &&
               (isEvaluatorResultTypeMappable<EvaluatorType>)) {
    constexpr auto has_data_method =
        requires(const EvaluatorType& rf, const cell_index<Space>& re,
                 const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, re, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::const_qt_ref<UnitType, real>(
          *(this->evaluator.data(unsafe, e, i)));
    } else {
      return tfel::math::const_qt_ref<UnitType, real>(
          *(this->evaluator(e, i).data()));
    }
  }

  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr auto QuantityModifier<EvaluatorType, UnitType>::operator()(
      const cell_workspace<Space>& wk,
      const cell_index<Space>& e,
      const quadrature_point_index<Space>& i) const
      requires((internals::EvaluatorResultQuery<EvaluatorType>::b4) &&
               (isEvaluatorResultTypeMappable<EvaluatorType>)) {
    constexpr auto has_data_method = requires(
        const EvaluatorType& rf, const cell_workspace<Space>& rwk,
        const cell_index<Space>& re, const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, rwk, re, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::const_qt_ref<UnitType, real>(
          *(this->evaluator.data(unsafe, wk, e, i)));
    } else {
      return tfel::math::const_qt_ref<UnitType, real>(
          *(this->evaluator(wk, e, i).data()));
    }
  }

  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  decltype(auto) getSpace(const QuantityModifier<EvaluatorType, UnitType>& e) {
    return e.getSpace();
  }  // end of getSpace

  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr bool check(AbstractErrorHandler& eh,
                       const QuantityModifier<EvaluatorType, UnitType>& e) {
    return e.check(eh);
  }  // end of check

  template <EvaluatorConcept EvaluatorType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr size_type getNumberOfComponents(
      const QuantityModifier<EvaluatorType, UnitType>& e) {
    return e.getNumberOfComponents();
  }  // end of getNumberOfComponents

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TFEL_QUANTITYMODIFIER!_IXX */
