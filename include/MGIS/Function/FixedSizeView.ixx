/*!
 * \file   FixedSizeView.ixx
 * \brief
 * \author th202608
 * \date   07/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_FIXEDSIZEVIEW_IXX
#define LIB_MGIS_FUNCTION_FIXEDSIZEVIEW_IXX

#include "MGIS/Raise.hxx"

namespace mgis::function {

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0) bool FixedSizeView<Space, N, is_mutable>::checkPreconditions(
      const FunctionView<Space, {}, is_mutable>& values) noexcept {
    return values.getNumberOfComponents() == N;
  }  // end of checkPreconditions

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0) FixedSizeView<Space, N, is_mutable>::FixedSizeView(
      const FunctionView<Space, {}, is_mutable>& values)
      : function(values) {
    raise_if(!checkPreconditions(values),
             "FixedSizeImmutableView::FixedSizeImmutableView: "
             "unmatched size");
  }  // end of FixedSizeView

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0) bool FixedSizeView<Space, N, is_mutable>::check(
      Context&) const noexcept {
    return true;
  }

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N >
           0) void FixedSizeView<Space, N, is_mutable>::allocateWorkspace() {}

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0) const Space& FixedSizeView<Space, N, is_mutable>::getSpace()
      const {
    return this->function.getSpace();
  }

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0) constexpr size_type
      FixedSizeView<Space, N, is_mutable>::getNumberOfComponents()
          const noexcept {
    return N;
  }

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0) auto FixedSizeView<Space, N, is_mutable>::operator()(
      const element_index<Space>& i) const
      requires(ElementSpaceConcept<Space> && !(hasElementWorkspace<Space>)) {
    if constexpr (N == 1) {
      return this->function.getValue(i);
    } else {
      return this->function.template getValues<N>(i);
    }
  }

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0) auto FixedSizeView<Space, N, is_mutable>::operator()(
      const element_workspace<Space>&, const element_index<Space>& i) const
      requires(ElementSpaceConcept<Space>&& hasElementWorkspace<Space>) {
    if constexpr (N == 1) {
      return this->function.getValue(i);
    } else {
      return this->function.template getValues<N>(i);
    }
  }

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0) auto FixedSizeView<Space, N, is_mutable>::operator()(
      const cell_index<Space> e, const quadrature_point_index<Space> i) const
      requires(QuadratureSpaceConcept<Space> && (!hasCellWorkspace<Space>)) {
    if constexpr (N == 1) {
      return this->function.getValue(e, i);
    } else {
      return this->function.template getValues<N>(e, i);
    }
  }

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0) auto FixedSizeView<Space, N, is_mutable>::operator()(
      const cell_workspace<Space>&,
      const cell_index<Space> e,
      const quadrature_point_index<Space> i) const
      requires(QuadratureSpaceConcept<Space>&& hasCellWorkspace<Space>) {
    if constexpr (N == 1) {
      return this->function.getValue(e, i);
    } else {
      return this->function.template getValues<N>(e, i);
    }
  }

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)
      typename FixedSizeView<Space, N, is_mutable>::mutable_value_type
      FixedSizeView<Space, N, is_mutable>::operator()(
          const element_index<Space>& i)  //
      requires(is_mutable&& ElementSpaceConcept<Space> &&
               !(hasElementWorkspace<Space>)) {
    if constexpr (N == 1) {
      return this->function.getValue(i);
    } else {
      return this->function.template getValues<N>(i);
    }
  }

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)
      typename FixedSizeView<Space, N, is_mutable>::mutable_value_type
      FixedSizeView<Space, N, is_mutable>::operator()(
          const element_workspace<Space>&,
          const element_index<Space>& i)  //
      requires(is_mutable&& ElementSpaceConcept<Space>&&
                   hasElementWorkspace<Space>) {
    if constexpr (N == 1) {
      return this->function.getValue(i);
    } else {
      return this->function.template getValues<N>(i);
    }
  }

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)
      typename FixedSizeView<Space, N, is_mutable>::mutable_value_type
      FixedSizeView<Space, N, is_mutable>::operator()(
          const cell_index<Space> e,
          const quadrature_point_index<Space> i)  //
      requires(is_mutable&& QuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    if constexpr (N == 1) {
      return this->function.getValue(e, i);
    } else {
      return this->function.template getValues<N>(e, i);
    }
  }

  template <FunctionalSpaceConcept Space, size_type N, bool is_mutable>
  requires(N > 0)
      typename FixedSizeView<Space, N, is_mutable>::mutable_value_type
      FixedSizeView<Space, N, is_mutable>::operator()(
          const cell_workspace<Space>&,
          const cell_index<Space> e,
          const quadrature_point_index<Space> i)  //
      requires(is_mutable&& QuadratureSpaceConcept<Space>&&
                   hasCellWorkspace<Space>) {
    if constexpr (N == 1) {
      return this->function.getValue(e, i);
    } else {
      return this->function.template getValues<N>(e, i);
    }
  }

  template <size_type N, FunctionalSpaceConcept Space, bool is_mutable>
  auto view(const FunctionView<Space, {}, is_mutable>& f) requires(N > 0) {
    return FixedSizeView<Space, N, false>(f);
  }  // end of view

  template <size_type N, FunctionalSpaceConcept Space>
  auto view(Function<Space, dynamic_extent>& f)  //
      requires((N > 0) && (N != dynamic_extent)) {
    return FixedSizeView<Space, N, true>(f.view());
  }  // end of view

  template <size_type N, FunctionalSpaceConcept Space>
  auto view(const Function<Space, dynamic_extent>& f)  //
      requires((N > 0) && (N != dynamic_extent)) {
    return FixedSizeView<Space, N, false>(f.view());
  }  // end of view

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_FIXEDSIZEVIEW_IXX */
