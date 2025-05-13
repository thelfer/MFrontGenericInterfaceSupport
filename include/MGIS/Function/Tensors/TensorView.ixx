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
    constexpr auto N = compile_time_size<TensorType>;
    return tfel::math::map<TensorType>(
        this->function.template getValues<N>(i).data());
  }

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  auto TensorView<Space, TensorType, is_mutable>::operator()(
      const element_workspace<Space>&, const element_index<Space>& i) const
      requires(ElementSpaceConcept<Space>&& hasElementWorkspace<Space>) {
    constexpr auto N = compile_time_size<TensorType>;
    return tfel::math::map<TensorType>(
        this->function.template getValues<N>(i).data());
  }

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  auto TensorView<Space, TensorType, is_mutable>::operator()(
      const cell_index<Space> e, const quadrature_point_index<Space> i) const
      requires(QuadratureSpaceConcept<Space> && (!hasCellWorkspace<Space>)) {
    constexpr auto N = compile_time_size<TensorType>;
    return tfel::math::map<TensorType>(
        this->function.template getValues<N>(e, i).data());
  }

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  auto TensorView<Space, TensorType, is_mutable>::operator()(
      const cell_workspace<Space>&,
      const cell_index<Space> e,
      const quadrature_point_index<Space> i) const
      requires(QuadratureSpaceConcept<Space>&& hasCellWorkspace<Space>) {
    constexpr auto N = compile_time_size<TensorType>;
    return tfel::math::map<TensorType>(
        this->function.template getValues<N>(e, i).data());
  }

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  auto TensorView<Space, TensorType, is_mutable>::operator()(
      const element_index<Space>& i)  //
      requires(is_mutable&& ElementSpaceConcept<Space> &&
               !(hasElementWorkspace<Space>)) {
    constexpr auto N = compile_time_size<TensorType>;
    return tfel::math::map<TensorType>(
        this->function.template getValues<N>(i).data());
  }

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  auto TensorView<Space, TensorType, is_mutable>::operator()(
      const element_workspace<Space>&,
      const element_index<Space>& i)  //
      requires(is_mutable&& ElementSpaceConcept<Space>&&
                   hasElementWorkspace<Space>) {
    constexpr auto N = compile_time_size<TensorType>;
    return tfel::math::map<TensorType>(
        this->function.template getValues<N>(i).data());
  }

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  auto TensorView<Space, TensorType, is_mutable>::operator()(
      const cell_index<Space> e,
      const quadrature_point_index<Space> i)  //
      requires(is_mutable&& QuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    constexpr auto N = compile_time_size<TensorType>;
    return tfel::math::map<TensorType>(
        this->function.template getValues<N>(e, i).data());
  }

  template <FunctionalSpaceConcept Space,
            TensorConcept TensorType,
            bool is_mutable>
  auto TensorView<Space, TensorType, is_mutable>::operator()(
      const cell_workspace<Space>&,
      const cell_index<Space> e,
      const quadrature_point_index<Space> i)  //
      requires(is_mutable&& QuadratureSpaceConcept<Space>&&
                   hasCellWorkspace<Space>) {
    constexpr auto N = compile_time_size<TensorType>;
    return tfel::math::map<TensorType>(
        this->function.template getValues<N>(e, i).data());
  }

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_TENSORVIEW_IXX */
