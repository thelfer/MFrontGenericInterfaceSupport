/*!
 * \file   MGIS/Function/TFEL/QuantityView.ixx
 * \brief
 * \author Thomas Helfer
 * \date   15/08/2026
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_TFEL_QUANTITYVIEW_IXX
#define LIB_MGIS_FUNCTION_TFEL_QUANTITYVIEW_IXX

namespace mgis::function {

  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr bool QuantityView<FunctionType, UnitType>::checkPreconditions(
      AbstractErrorHandler& eh, const FunctionType& values) {
    if (internals::disambiguateGetNumberOfComponents(values) != 1) {
      return eh.registerErrorMessage("invalid number of components");
    }
    return true;
  }  // end of checkPreconditions

  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr QuantityView<FunctionType, UnitType>::QuantityView(
      FunctionType& values)
      : QuantityView(preconditions_check, values) {}  // end of QuantityView

  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  template <bool doPreconditionsCheck>
  constexpr QuantityView<FunctionType, UnitType>::QuantityView(
      const PreconditionsCheck<doPreconditionsCheck>& pcheck,
      FunctionType& values)
      : PreconditionsChecker<QuantityView>(pcheck, values),
        function(make_view(values)) {}  // end of QuantityView

  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr bool QuantityView<FunctionType, UnitType>::check(
      AbstractErrorHandler& eh) const {
    if (internals::disambiguateGetNumberOfComponents(this->function) != 1) {
      return eh.registerErrorMessage("invalid number of components");
    }
    return true;
  }  // end of check

  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr decltype(auto) QuantityView<FunctionType, UnitType>::getSpace()
      const {
    return internals::disambiguateGetSpace(this->function);
  }

  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr size_type
  QuantityView<FunctionType, UnitType>::getNumberOfComponents() const noexcept {
    return 1;
  }  // end of getNumberOfComponents

  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr auto QuantityView<FunctionType, UnitType>::operator()(
      const element_index<Space>& i) const
      requires((internals::FunctionResultQuery<FunctionType>::b1) &&
               (isFunctionConstResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method =
        requires(const FunctionType& rf, const element_index<Space>& ri) {
      { rf.data(unsafe, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return ::tfel::math::const_qt_ref<UnitType, real>(
          *(this->function.data(unsafe, i)));
    } else {
      return ::tfel::math::const_qt_ref<UnitType, real>(
          *(this->function(i).data()));
    }
  }

  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr auto QuantityView<FunctionType, UnitType>::operator()(
      const element_workspace<Space>& wk, const element_index<Space>& i) const
      requires((internals::FunctionResultQuery<FunctionType>::b2) &&
               (isFunctionConstResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method =
        requires(const FunctionType& rf, const element_workspace<Space>& rwk,
                 const element_index<Space>& ri) {
      { rf.data(unsafe, rwk, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return ::tfel::math::const_qt_ref<UnitType, real>(
          *(this->function.data(unsafe, wk, i)));
    } else {
      return ::tfel::math::const_qt_ref<UnitType, real>(
          *(this->function(wk, i).data()));
    }
  }

  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr auto QuantityView<FunctionType, UnitType>::operator()(
      const cell_index<Space>& e, const quadrature_point_index<Space>& i) const
      requires((internals::FunctionResultQuery<FunctionType>::b3) &&
               (isFunctionConstResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method =
        requires(const FunctionType& rf, const cell_index<Space>& re,
                 const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, re, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return ::tfel::math::const_qt_ref<UnitType, real>(
          *(this->function.data(unsafe, e, i)));
    } else {
      return ::tfel::math::const_qt_ref<UnitType, real>(
          *(this->function(e, i).data()));
    }
  }

  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr auto QuantityView<FunctionType, UnitType>::operator()(
      const cell_workspace<Space>& wk,
      const cell_index<Space>& e,
      const quadrature_point_index<Space>& i) const
      requires((internals::FunctionResultQuery<FunctionType>::b4) &&
               (isFunctionConstResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method = requires(
        const FunctionType& rf, const cell_workspace<Space>& rwk,
        const cell_index<Space>& re, const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, rwk, re, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return ::tfel::math::const_qt_ref<UnitType, real>(
          *(this->function.data(unsafe, wk, e, i)));
    } else {
      return ::tfel::math::const_qt_ref<UnitType, real>(
          *(this->function(wk, e, i).data()));
    }
  }

  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr auto QuantityView<FunctionType, UnitType>::operator()(
      const element_index<Space>& i)  //
      requires((internals::FunctionResultQuery<FunctionType>::b1) &&
               (isFunctionResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method =
        requires(FunctionType & rf, const element_index<Space>& ri) {
      { rf.data(unsafe, ri) } -> std::same_as<real*>;
    };
    if constexpr (has_data_method) {
      return ::tfel::math::qt_ref<UnitType, real>(
          *(this->function.data(unsafe, i)));
    } else {
      return ::tfel::math::qt_ref<UnitType, real>(*(this->function(i).data()));
    }
  }

  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr auto QuantityView<FunctionType, UnitType>::operator()(
      const element_workspace<Space>& wk,
      const element_index<Space>& i)  //
      requires((internals::FunctionResultQuery<FunctionType>::b2) &&
               (isFunctionResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method =
        requires(FunctionType & rf, const element_workspace<Space>& rwk,
                 const element_index<Space>& ri) {
      { rf.data(unsafe, rwk, ri) } -> std::same_as<real*>;
    };
    if constexpr (has_data_method) {
      return ::tfel::math::qt_ref<UnitType, real>(
          *(this->function.data(unsafe, wk, i)));
    } else {
      return ::tfel::math::qt_ref<UnitType, real>(
          *(this->function(wk, i).data()));
    }
  }

  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr auto QuantityView<FunctionType, UnitType>::operator()(
      const cell_index<Space>& e,
      const quadrature_point_index<Space>& i)  //
      requires((internals::FunctionResultQuery<FunctionType>::b3) &&
               (isFunctionResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method =
        requires(FunctionType & rf, const cell_index<Space>& re,
                 const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, re, ri) } -> std::same_as<real*>;
    };
    if constexpr (has_data_method) {
      return ::tfel::math::qt_ref<UnitType, real>(
          *(this->function.data(unsafe, e, i)));
    } else {
      return ::tfel::math::qt_ref<UnitType, real>(
          *(this->function(e, i).data()));
    }
  }

  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr auto QuantityView<FunctionType, UnitType>::operator()(
      const cell_workspace<Space>& wk,
      const cell_index<Space>& e,
      const quadrature_point_index<Space>& i)  //
      requires((internals::FunctionResultQuery<FunctionType>::b4) &&
               (isFunctionResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method = requires(
        FunctionType & rf, const cell_workspace<Space>& rwk,
        const cell_index<Space>& re, const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, rwk, re, ri) } -> std::same_as<real*>;
    };
    if constexpr (has_data_method) {
      return ::tfel::math::qt_ref<UnitType, real>(
          *(this->function.data(unsafe, wk, e, i)));
    } else {
      return ::tfel::math::qt_ref<UnitType, real>(
          *(this->function(wk, e, i).data()));
    }
  }

  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr decltype(auto) getSpace(
      const QuantityView<FunctionType, UnitType>& v) {
    return v.getSpace();
  }  // end of getSpace

  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr bool check(AbstractErrorHandler& eh,
                       const QuantityView<FunctionType, UnitType>& v) {
    return v.check(eh);
  }  // end of check

  template <FunctionConcept FunctionType,
            ::tfel::math::unit::UnitConcept UnitType>
  constexpr mgis::size_type getNumberOfComponents(
      const QuantityView<FunctionType, UnitType>& v) noexcept {
    return v.getNumberOfComponents();
  }  // end of getNumberOfComponents

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TFEL_QUANTITYVIEW_IXX */
