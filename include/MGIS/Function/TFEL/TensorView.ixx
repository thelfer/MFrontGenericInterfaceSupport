
/*!
 * \file   MGIS/Function/TFEL/TensorView.ixx
 * \brief
 * \author Thomas Helfer
 * \date   10/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_TFEL_TENSORVIEW_IXX
#define LIB_MGIS_FUNCTION_TFEL_TENSORVIEW_IXX

#include <utility>

namespace mgis::function {

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr bool TensorView<FunctionType, TensorType>::checkPreconditions(
      AbstractErrorHandler& eh, const FunctionType& values) {
    const auto nc = internals::disambiguateGetNumberOfComponents(values);
    if (nc != compile_time_size<TensorType>) {
      return eh.registerErrorMessage("invalid number of components");
    }
    return true;
  }  // end of checkPreconditions

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr TensorView<FunctionType, TensorType>::TensorView(
      FunctionType& values)
      : TensorView(preconditions_check, values) {}  // end of TensorView

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  template <bool doPreconditionsCheck>
  constexpr TensorView<FunctionType, TensorType>::TensorView(
      const PreconditionsCheck<doPreconditionsCheck>& pcheck,
      FunctionType& values)
      : PreconditionsChecker<TensorView>(pcheck, values),
        function(make_view(values)) {}  // end of TensorView

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr bool TensorView<FunctionType, TensorType>::check(
      AbstractErrorHandler& eh) const {
    const auto nc =
        internals::disambiguateGetNumberOfComponents(this->function);
    if (nc != compile_time_size<TensorType>) {
      return eh.registerErrorMessage("invalid number of components");
    }
    return true;
  }  // end of check

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr const typename TensorView<FunctionType, TensorType>::Space&
  TensorView<FunctionType, TensorType>::getSpace() const {
    return internals::disambiguateGetSpace(this->function);
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr size_type
  TensorView<FunctionType, TensorType>::getNumberOfComponents() const noexcept {
    return compile_time_size<TensorType>;
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr auto TensorView<FunctionType, TensorType>::operator()(
      const element_index<Space>& i) const
      requires((internals::FunctionResultQuery<FunctionType>::b1) &&
               (isFunctionConstResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method =
        requires(const FunctionType& rf, const element_index<Space>& ri) {
      { rf.data(unsafe, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<const TensorType>(this->function.data(unsafe, i));
    } else {
      return tfel::math::map<const TensorType>(this->function(i).data());
    }
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr auto TensorView<FunctionType, TensorType>::operator()(
      const element_workspace<Space>& wk, const element_index<Space>& i) const
      requires((internals::FunctionResultQuery<FunctionType>::b2) &&
               (isFunctionConstResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method =
        requires(const FunctionType& rf, const element_workspace<Space>& rwk,
                 const element_index<Space>& ri) {
      { rf.data(unsafe, rwk, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<const TensorType>(
          this->function.data(unsafe, wk, i));
    } else {
      return tfel::math::map<const TensorType>(this->function(wk, i).data());
    }
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr auto TensorView<FunctionType, TensorType>::operator()(
      const cell_index<Space>& e, const quadrature_point_index<Space>& i) const
      requires((internals::FunctionResultQuery<FunctionType>::b3) &&
               (isFunctionConstResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method =
        requires(const FunctionType& rf, const cell_index<Space> re,
                 const quadrature_point_index<Space> ri) {
      { rf.data(unsafe, re, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<const TensorType>(
          this->function.data(unsafe, e, i));
    } else {
      return tfel::math::map<const TensorType>(this->function(e, i).data());
    }
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr auto TensorView<FunctionType, TensorType>::operator()(
      const cell_workspace<Space>& wk,
      const cell_index<Space>& e,
      const quadrature_point_index<Space>& i) const
      requires((internals::FunctionResultQuery<FunctionType>::b4) &&
               (isFunctionConstResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method = requires(
        const FunctionType& rf, const cell_workspace<Space>& rwk,
        const cell_index<Space> re, const quadrature_point_index<Space> ri) {
      { rf.data(unsafe, rwk, re, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<const TensorType>(
          this->function.data(unsafe, wk, e, i));
    } else {
      return tfel::math::map<const TensorType>(this->function(wk, e, i).data());
    }
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr auto TensorView<FunctionType, TensorType>::operator()(
      const element_index<Space>& i)  //
      requires((internals::FunctionResultQuery<FunctionType>::b1) &&
               (isFunctionResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method =
        requires(FunctionType & rf, const element_index<Space>& ri) {
      { rf.data(unsafe, ri) } -> std::same_as<real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<TensorType>(this->function.data(unsafe, i));
    } else {
      return tfel::math::map<TensorType>(this->function(i).data());
    }
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr auto TensorView<FunctionType, TensorType>::operator()(
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
      return tfel::math::map<TensorType>(this->function.data(unsafe, wk, i));
    } else {
      return tfel::math::map<TensorType>(this->function(wk, i).data());
    }
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr auto TensorView<FunctionType, TensorType>::operator()(
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
      return tfel::math::map<TensorType>(this->function.data(unsafe, e, i));
    } else {
      return tfel::math::map<TensorType>(this->function(e, i).data());
    }
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr auto TensorView<FunctionType, TensorType>::operator()(
      const cell_workspace<Space>& wk,
      const cell_index<Space>& e,
      const quadrature_point_index<Space>& i)  //
      requires((internals::FunctionResultQuery<FunctionType>::b4) &&
               (isFunctionResultTypeMappable<FunctionType>)) {
    constexpr auto has_data_method = requires(
        FunctionType & rf, const cell_workspace<Space>& rwk,
        const cell_index<Space>& re, const quadrature_point_index<Space>& ri) {
      { rf.data(unsafe, rwk, ri) } -> std::same_as<real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<TensorType>(this->function.data(unsafe, wk, e, i));
    } else {
      return tfel::math::map<TensorType>(this->function(wk, e, i).data());
    }
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr decltype(auto) getSpace(
      const TensorView<FunctionType, TensorType>& t) {
    return t.getSpace();
  }

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr bool check(AbstractErrorHandler& eh,
                       const TensorView<FunctionType, TensorType>& v) {
    return v.check(eh);
  }  // end of check

  template <FunctionConcept FunctionType, TensorConcept TensorType>
  constexpr size_type getNumberOfComponents(
      const TensorView<FunctionType, TensorType>& v) noexcept {
    return v.getNumberOfComponents();
  }  // end of getNumberOfComponents

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TFEL_TENSORVIEW_IXX */
