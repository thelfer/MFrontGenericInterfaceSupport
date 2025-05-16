/*!
 * \file   MGIS/Function/TensorView.ixx
 * \brief
 * \author Thomas Helfer
 * \date   10/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_TENSORVIEW_IXX
#define LIB_MGIS_FUNCTION_TENSORVIEW_IXX

#include <utility>

namespace mgis::function {

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  bool TensorView<Space, TensorType, is_mutable>::checkPreconditions(
      const FunctionView<Space, {}, is_mutable>& values) noexcept {
    return values.getNumberOfComponents() == compile_time_size<TensorType>;
  }  // end of checkPreconditions

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  TensorView<Space, TensorType, is_mutable>::TensorView(
      const FunctionView<Space, {}, is_mutable>& values)
      : function(values) {
    raise_if(!checkPreconditions(values),
             "FixedSizeImmutableView::FixedSizeImmutableView: "
             "unmatched size");
  }  // end of FixedSizeView

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  bool TensorView<Space, TensorType, is_mutable>::check(
      Context&) const noexcept {
    return checkPreconditions(this->function);
  }

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  void TensorView<Space, TensorType, is_mutable>::allocateWorkspace() {}

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  const Space& TensorView<Space, TensorType, is_mutable>::getSpace() const {
    return this->function.getSpace();
  }

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  constexpr size_type
  TensorView<Space, TensorType, is_mutable>::getNumberOfComponents()
      const noexcept {
    return compile_time_size<TensorType>;
  }

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  auto TensorView<Space, TensorType, is_mutable>::operator()(
      const element_index<Space>& i) const
      requires(ElementSpaceConcept<Space> && !(hasElementWorkspace<Space>)) {
    constexpr auto has_data_method =
        requires(const FunctionView<Space, {}, is_mutable>& rf,
                 const element_index<Space>& ri) {
      { rf.data(ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<TensorType>(this->function.data(unsafe, i));
    } else {
      return tfel::math::map<TensorType>(this->function(i).data());
    }
  }

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  auto TensorView<Space, TensorType, is_mutable>::operator()(
      const element_workspace<Space>& wk, const element_index<Space>& i) const
      requires(ElementSpaceConcept<Space>&& hasElementWorkspace<Space>) {
    constexpr auto has_data_method = requires(
        const FunctionView<Space, {}, is_mutable>& rf,
        const element_workspace<Space>& rwk, const element_index<Space>& ri) {
      { rf.data(rwk, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<TensorType>(this->function.data(unsafe, wk, i));
    } else {
      return tfel::math::map<TensorType>(this->function(wk, i).data());
    }
  }

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  auto TensorView<Space, TensorType, is_mutable>::operator()(
      const cell_index<Space>& e, const quadrature_point_index<Space>& i) const
      requires(QuadratureSpaceConcept<Space> && (!hasCellWorkspace<Space>)) {
    constexpr auto has_data_method = requires(
        const FunctionView<Space, {}, is_mutable>& rf,
        const cell_index<Space> re, const quadrature_point_index<Space> ri) {
      { rf.data(re, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<TensorType>(this->function.data(unsafe, e, i));
    } else {
      return tfel::math::map<TensorType>(this->function(e, i).data());
    }
  }

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  auto TensorView<Space, TensorType, is_mutable>::operator()(
      const cell_workspace<Space>& wk,
      const cell_index<Space>& e,
      const quadrature_point_index<Space>& i) const
      requires(QuadratureSpaceConcept<Space>&& hasCellWorkspace<Space>) {
    constexpr auto has_data_method =
        requires(const FunctionView<Space, {}, is_mutable>& rf,
                 const cell_workspace<Space>& rwk, const cell_index<Space> re,
                 const quadrature_point_index<Space> ri) {
      { rf.data(rwk, re, ri) } -> std::same_as<const real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<TensorType>(this->function.data(unsafe, wk, e, i));
    } else {
      return tfel::math::map<TensorType>(this->function(wk, e, i).data());
    }
  }

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  auto TensorView<Space, TensorType, is_mutable>::operator()(
      const element_index<Space>& i)  //
      requires(is_mutable&& ElementSpaceConcept<Space> &&
               !(hasElementWorkspace<Space>)) {
    constexpr auto has_data_method =
        requires(FunctionView<Space, {}, is_mutable> & rf,
                 const element_index<Space>& ri) {
      { rf.data(ri) } -> std::same_as<real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<TensorType>(this->function.data(unsafe, i));
    } else {
      return tfel::math::map<TensorType>(this->function(i).data());
    }
  }

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  auto TensorView<Space, TensorType, is_mutable>::operator()(
      const element_workspace<Space>& wk,
      const element_index<Space>& i)  //
      requires(is_mutable&& ElementSpaceConcept<Space>&&
                   hasElementWorkspace<Space>) {
    constexpr auto has_data_method = requires(
        FunctionView<Space, {}, is_mutable> & rf,
        const element_workspace<Space>& rwk, const element_index<Space>& ri) {
      { rf.data(rwk, ri) } -> std::same_as<real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<TensorType>(this->function.data(unsafe, wk, i));
    } else {
      return tfel::math::map<TensorType>(this->function(wk, i).data());
    }
  }

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  auto TensorView<Space, TensorType, is_mutable>::operator()(
      const cell_index<Space>& e,
      const quadrature_point_index<Space>& i)  //
      requires(is_mutable&& QuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    constexpr auto has_data_method = requires(
        FunctionView<Space, {}, is_mutable> & rf, const cell_index<Space>& re,
        const quadrature_point_index<Space>& ri) {
      { rf.data(re, ri) } -> std::same_as<real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<TensorType>(this->function.data(unsafe, e, i));
    } else {
      return tfel::math::map<TensorType>(this->function(e, i).data());
    }
  }

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  auto TensorView<Space, TensorType, is_mutable>::operator()(
      const cell_workspace<Space>& wk,
      const cell_index<Space>& e,
      const quadrature_point_index<Space>& i)  //
      requires(is_mutable&& QuadratureSpaceConcept<Space>&&
                   hasCellWorkspace<Space>) {
    constexpr auto has_data_method =
        requires(FunctionView<Space, {}, is_mutable> & rf,
                 const cell_workspace<Space>& rwk, const cell_index<Space>& re,
                 const quadrature_point_index<Space>& ri) {
      { rf.data(rwk, ri) } -> std::same_as<real*>;
    };
    if constexpr (has_data_method) {
      return tfel::math::map<TensorType>(this->function.data(unsafe, wk, i));
    } else {
      return tfel::math::map<TensorType>(this->function(wk, i).data());
    }
  }

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TENSORVIEW_IXX */
