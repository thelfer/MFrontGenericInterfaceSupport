/*!
 * \file   MGIS/Function/Function.ixx
 * \brief
 * \author Thomas Helfer
 * \date   10/03/2021
 */

#ifndef LIB_MGIS_FUNCTION_FUNCTION_IXX
#define LIB_MGIS_FUNCTION_FUNCTION_IXX

#include <concepts>
#include "MGIS/Raise.hxx"

namespace mgis::function {

  template <size_type data_size>
  constexpr bool FunctionDataSize<data_size>::isScalar() const noexcept {
    return data_size == 1;
  }  // end of isScalar

  template <size_type data_size>
  constexpr size_type FunctionDataSize<data_size>::getNumberOfComponents()
      const noexcept {
    return data_size;
  }  // end of FunctionDataSize<dynamic_extent>::getNumberOfComponents

  inline bool FunctionDataSize<dynamic_extent>::isScalar() const noexcept {
    return this->getNumberOfComponents() == 1;
  }  // end of isScalar

  inline size_type FunctionDataSize<dynamic_extent>::getNumberOfComponents()
      const noexcept {
    return this->data_size;
  }  // end of FunctionDataSize<dynamic_extent>::getNumberOfComponents

  template <size_type data_stride>
  constexpr size_type FunctionDataStride<data_stride>::getDataStride()
      const noexcept {
    return data_stride;
  }  // end of getDataStride

  inline size_type FunctionDataStride<dynamic_extent>::getDataStride()
      const noexcept {
    return this->data_stride;
  }  // end of getDataStride

  constexpr bool has_dynamic_properties(const DataLayoutDescription& layout) {
    return (layout.size == dynamic_extent) ||  //
           (layout.stride == dynamic_extent);
  }  // end of has_dynamic_properties

  template <DataLayoutDescription layout>
  inline size_type DataLayout<layout>::getDataOffset(
      const size_type i) const noexcept {
    if constexpr (layout.stride != dynamic_extent) {
      return i * layout.stride;
    } else {
      return i * (this->getDataStride());
    }
  }  // end of getDataOffset

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  bool FunctionView<Space, layout, is_mutable>::checkPreconditions(
      std::shared_ptr<const Space> s,
      typename FunctionView::ExternalData v) requires((layout.size !=
                                                       dynamic_extent) &&
                                                      (layout.stride !=
                                                       dynamic_extent)) {
    const auto space_size = s->size();
    if (space_size == 0) {
      // this may happen due to partionning in parallel
      return true;
    }
    return v.size() >= layout.stride * (space_size - 1) + layout.size;
  }

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  FunctionView<Space, layout, is_mutable>::FunctionView(
      std::shared_ptr<const Space> s,
      typename FunctionView::ExternalData v) requires((layout.size !=
                                                       dynamic_extent) &&
                                                      (layout.stride !=
                                                       dynamic_extent))
      : space(s), values(v) {
    raise_if(!checkPreconditions(s, v),
             "FunctionView::FunctionView:"
             " invalid values size");
  }

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  bool FunctionView<Space, layout, is_mutable>::checkPreconditions(
      std::shared_ptr<const Space> s,
      FunctionView::ExternalData v,
      const size_type dsize) requires((layout.size == dynamic_extent) &&
                                      (layout.stride != dynamic_extent)) {
    if (dsize <= 0) {
      return false;
    }
    const auto space_size = s->size();
    if (space_size == 0) {
      // this may happen due to partionning in parallel
      return true;
    }
    const auto min_size = (layout.stride) * (space_size - 1) + dsize;
    return v.size() >= min_size;
  }  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  FunctionView<Space, layout, is_mutable>::FunctionView(
      std::shared_ptr<const Space> s,
      FunctionView::ExternalData v,
      const size_type dsize) requires((layout.size == dynamic_extent) &&
                                      (layout.stride != dynamic_extent))
      : space(s), values(v) {
    raise_if(!checkPreconditions(s, v, dsize),
             "FunctionView::FunctionView:"
             " invalid values size");
    this->data_size = dsize;
  }  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  bool FunctionView<Space, layout, is_mutable>::checkPreconditions(
      std::shared_ptr<const Space> s,
      FunctionView::ExternalData v,
      const size_type dstride) requires((layout.size != dynamic_extent) &&
                                        (layout.stride == dynamic_extent)) {
    if (dstride <= 0) {
      return false;
    }
    const auto space_size = s->size();
    if (space_size == 0) {
      // this may happen due to partionning in parallel
      return true;
    }
    const auto min_size = dstride * (space_size - 1) + layout.size;
    return v.size() >= min_size;
  }  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  FunctionView<Space, layout, is_mutable>::FunctionView(
      std::shared_ptr<const Space> s,
      FunctionView::ExternalData v,
      const size_type dstride) requires((layout.size != dynamic_extent) &&
                                        (layout.stride == dynamic_extent))
      : space(s), values(v) {
    raise_if(!checkPreconditions(s, v, dstride),
             "FunctionView::FunctionView:"
             " invalid values size");
    this->data_stride = dstride;
  }  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  bool FunctionView<Space, layout, is_mutable>::checkPreconditions(
      std::shared_ptr<const Space> s,
      FunctionView::ExternalData v,
      const size_type dsize,
      const size_type dstride) requires((layout.size == dynamic_extent) &&
                                        (layout.stride == dynamic_extent)) {
    if (dsize <= 0) {
      return false;
    }
    if (dstride <= 0) {
      return false;
    }
    if (dsize > dstride) {
      return false;
    }
    const auto space_size = s->size();
    if (space_size == 0) {
      // this may happen due to partionning in parallel
      return true;
    }
    const auto min_size = dstride * (space_size - 1) + dsize;
    return v.size() >= min_size;
  }  // end of checkPreconditions

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  bool FunctionView<Space, layout, is_mutable>::checkPreconditions(
      std::shared_ptr<const Space> s,
      FunctionView::ExternalData v,
      const size_type dsize) requires((layout.size == dynamic_extent) &&
                                      (layout.stride == dynamic_extent)) {
    return checkPreconditions(s, v, dsize, dsize);
  }  // end of checkPreconditions

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  FunctionView<Space, layout, is_mutable>::FunctionView(
      std::shared_ptr<const Space> s,
      FunctionView::ExternalData v,
      const size_type dsize,
      const size_type dstride) requires((layout.size == dynamic_extent) &&
                                        (layout.stride == dynamic_extent))
      : space(s), values(v) {
    raise_if(!checkPreconditions(s, v, dsize, dstride),
             "FunctionView::FunctionView:"
             " invalid values size");
    this->data_size = dsize;
    this->data_stride = dstride;
  }  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  FunctionView<Space, layout, is_mutable>::FunctionView(
      std::shared_ptr<const Space> s,
      FunctionView::ExternalData v,
      const size_type dsize) requires((layout.size == dynamic_extent) &&
                                      (layout.stride == dynamic_extent))
      : FunctionView(s, v, dsize, dsize) {}  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  FunctionView<Space, layout, is_mutable>::FunctionView(
      std::shared_ptr<const Space> s,
      FunctionView::ExternalData v,
      const DataLayout<layout>& l) requires((layout.size == dynamic_extent) &&
                                            (layout.stride == dynamic_extent))
      : FunctionView(s, v, l.size, l.stride) {}  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  FunctionView<Space, layout, is_mutable>::FunctionView(
      std::shared_ptr<const Space> s,
      FunctionView::ExternalData v,
      const DataLayout<layout>& l) requires((layout.size != dynamic_extent) &&
                                            (layout.stride == dynamic_extent))
      : FunctionView(s, v, l.stride) {}

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  FunctionView<Space, layout, is_mutable>::FunctionView(
      std::shared_ptr<const Space> s,
      FunctionView::ExternalData v,
      const DataLayout<layout>& l) requires((layout.size == dynamic_extent) &&
                                            (layout.stride != dynamic_extent))
      : FunctionView(s, v, l.size) {}  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  bool FunctionView<Space, layout, is_mutable>::checkCompatibility(
      const FunctionView& v) const {
    if (this->getSpacePointer() != v.getSpacePointer()) {
      return false;
    }
    return this->data_size == v.getNumberOfComponents();
  }  // end of checkCompatibility

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  const Space& FunctionView<Space, layout, is_mutable>::getSpace()
      const noexcept {
    return *(this->space);
  }  // end of getSpace

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  std::shared_ptr<const Space>
  FunctionView<Space, layout, is_mutable>::getSpacePointer() const noexcept {
    return this->space;
  }  // end of getSpacePointer

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  bool FunctionView<Space, layout, is_mutable>::check(Context&) const noexcept {
    return true;
  }  // end of allocateWorkspace

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  void FunctionView<Space, layout, is_mutable>::allocateWorkspace() noexcept {
  }  // end of allocateWorkspace

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  real&
  FunctionView<Space, layout, is_mutable>::getValue(const size_type o) requires(
      allowScalarAccessor&& is_mutable&& LinearElementSpaceConcept<Space> &&
      (!hasElementWorkspace<Space>)) {
    return *(this->values.data() + this->getDataOffset(o));
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  typename FunctionView<Space, layout, is_mutable>::ValuesView
  FunctionView<Space, layout, is_mutable>::getValues(
      const size_type
          o) requires(is_mutable&& LinearElementSpaceConcept<Space> &&
                      (!hasElementWorkspace<Space>)) {
    return ValuesView(this->values.data() + this->getDataOffset(o),
                      this->getNumberOfComponents());
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  template <size_type N>
  std::span<real, N> FunctionView<Space, layout, is_mutable>::getValues(
      const size_type o) requires((layout.size == dynamic_extent) &&
                                  is_mutable &&
                                  LinearElementSpaceConcept<Space> &&
                                  (!hasElementWorkspace<Space>)) {
    return std::span<real, N>(this->values.data() + this->getDataOffset(o),
                              this->getNumberOfComponents());
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  real& FunctionView<Space, layout, is_mutable>::getValue(
      const size_type e,
      const size_type i) requires(allowScalarAccessor&& is_mutable&&
                                      LinearQuadratureSpaceConcept<Space> &&
                                  (!hasCellWorkspace<Space>)) {
    return this->getValue(this->space->getQuadraturePointOffset(e, i));
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  typename FunctionView<Space, layout, is_mutable>::ValuesView
  FunctionView<Space, layout, is_mutable>::getValues(
      const size_type e,
      const size_type
          i) requires(is_mutable&& LinearQuadratureSpaceConcept<Space> &&
                      (!hasCellWorkspace<Space>)) {
    return this->getValues(this->space->getQuadraturePointOffset(e, i));
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  template <size_type N>
  std::span<real, N> FunctionView<Space, layout, is_mutable>::getValues(
      const size_type e,
      const size_type i) requires((layout.size == dynamic_extent) &&
                                  is_mutable &&
                                  LinearQuadratureSpaceConcept<Space> &&
                                  (!hasCellWorkspace<Space>)) {
    return this->template getValues<N>(
        this->space->getQuadraturePointOffset(e, i));
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  const real& FunctionView<Space, layout, is_mutable>::getValue(
      const size_type o) const
      requires(allowScalarAccessor&& LinearElementSpaceConcept<Space> &&
               (!hasElementWorkspace<Space>)) {
    return *(this->values.data() + this->getDataOffset(o));
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  typename FunctionView<Space, layout, is_mutable>::ConstValuesView
  FunctionView<Space, layout, is_mutable>::getValues(const size_type o) const
      requires(LinearElementSpaceConcept<Space> &&
      (!hasElementWorkspace<Space>)) {
    return ConstValuesView(this->values.data() + this->getDataOffset(o),
                           this->getNumberOfComponents());
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  template <size_type N>
  std::span<const real, N> FunctionView<Space, layout, is_mutable>::getValues(
      const size_type o) const requires(LinearElementSpaceConcept<Space> &&
      (!hasElementWorkspace<Space>)) {
    return std::span<const real, N>(
        this->values.data() + this->getDataOffset(o),
        this->getNumberOfComponents());
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  const real& FunctionView<Space, layout, is_mutable>::getValue(
      const size_type e, const size_type i) const
      requires(allowScalarAccessor&& LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return this->getValue(this->space->getQuadraturePointOffset(e, i));
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  typename FunctionView<Space, layout, is_mutable>::ConstValuesView
  FunctionView<Space, layout, is_mutable>::getValues(const size_type e,
                                                     const size_type i) const
      requires(LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return this->getValues(this->space->getQuadraturePointOffset(e, i));
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  template <size_type N>
  std::span<const real, N> FunctionView<Space, layout, is_mutable>::getValues(
      const size_type e, const size_type i) const
      requires(LinearQuadratureSpaceConcept<Space> &&
      (!hasCellWorkspace<Space>)) {
    return this->template getValues<N>(
        this->space->getQuadraturePointOffset(e, i));
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  typename FunctionView<Space, layout, is_mutable>::ValuesView
  FunctionView<Space, layout, is_mutable>::operator()(const size_type o)  //
      requires(is_mutable&& LinearElementSpaceConcept<Space> &&
               (!hasElementWorkspace<Space>)) {
    return ValuesView(this->values.data() + this->getDataOffset(o),
                           this->getNumberOfComponents());
  }  // end of operator()

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  typename FunctionView<Space, layout, is_mutable>::ValuesView
  FunctionView<Space, layout, is_mutable>::operator()(const size_type e,
                                                      const size_type i)  //
      requires(is_mutable&& LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return this->operator()(this->space->getQuadraturePointOffset(e, i));
  }  // end of operator()

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  typename FunctionView<Space, layout, is_mutable>::ConstValuesView
  FunctionView<Space, layout, is_mutable>::operator()(const size_type o) const
      requires(LinearElementSpaceConcept<Space> &&
               (!hasElementWorkspace<Space>)) {
    return ConstValuesView(this->values.data() + this->getDataOffset(o),
                           this->getNumberOfComponents());
  }  // end of operator()

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  typename FunctionView<Space, layout, is_mutable>::ConstValuesView
  FunctionView<Space, layout, is_mutable>::operator()(const size_type e,
                                                      const size_type i) const
      requires(LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return this->operator()(this->space->getQuadraturePointOffset(e, i));
  }  // end of operator()

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  std::span<const real> FunctionView<Space, layout, is_mutable>::getValues()
      const {
    return this->values;
  }

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  FunctionView<Space, layout, is_mutable>::~FunctionView() noexcept = default;

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_FUNCTION_IXX */
