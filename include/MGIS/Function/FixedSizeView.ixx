/*!
 * \file   FixedSizeView.ixx
 * \brief
 * \author th202608
 * \date   07/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_FIXEDSIZEEVALUATOR_IXX
#define LIB_MGIS_FUNCTION_FIXEDSIZEEVALUATOR_IXX

#include "MGIS/Raise.hxx"

namespace mgis::function {

  template <FunctionalSpaceConcept Space, size_type N>
  requires(N > 0) bool FixedSizeView<Space, N>::checkPreconditions(
      const ImmutableFunctionView<Space, {}>& values) noexcept {
    return values.getNumberOfComponents() == N;
  }  // end of checkPreconditions

  template <FunctionalSpaceConcept Space, size_type N>
  requires(N > 0) FixedSizeView<Space, N>::FixedSizeView(
      const ImmutableFunctionView<Space, {}>& values)
      : function(values) {
    raise_if(!checkPreconditions(values),
             "FixedSizeImmutableView::FixedSizeImmutableView: "
             "unmatched size");
  }  // end of FixedSizeView

  template <FunctionalSpaceConcept Space, size_type N>
  requires(N > 0) bool FixedSizeView<Space, N>::check(Context&) const noexcept {
    return true;
  }

  template <FunctionalSpaceConcept Space, size_type N>
  requires(N > 0) void FixedSizeView<Space, N>::allocateWorkspace() {}

  template <FunctionalSpaceConcept Space, size_type N>
  requires(N > 0) const Space& FixedSizeView<Space, N>::getSpace() const {
    return this->function.getSpace();
  }

  template <FunctionalSpaceConcept Space, size_type N>
  requires(N > 0) constexpr size_type
      FixedSizeView<Space, N>::getNumberOfComponents() const noexcept {
    return N;
  }

  template <FunctionalSpaceConcept Space, size_type N>
  requires(N > 0) auto FixedSizeView<Space, N>::operator()(
      const element_index<Space>& i) const
      requires(ElementSpaceConcept<Space> && !(hasElementWorkspace<Space>)) {
    if constexpr (N == 1) {
      return this->function.getValue(i);
    } else {
      return this->function.template getValues<N>(i);
    }
  }

  template <FunctionalSpaceConcept Space, size_type N>
  requires(N > 0) auto FixedSizeView<Space, N>::operator()(
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

  template <size_type N, FunctionalSpaceConcept Space>
  auto view(const ImmutableFunctionView<Space, {}>& f) requires(N > 0) {
    return FixedSizeView<Space, N>(f);
  }  // end of view

  template <size_type N, FunctionalSpaceConcept Space>
  auto view(const Function<Space, dynamic_extent>& f)  //
      requires((N > 0) && (N != dynamic_extent)) {
    return FixedSizeView<Space, N>(f.view());
  }  // end of view

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_FIXEDSIZEEVALUATOR_IXX */
