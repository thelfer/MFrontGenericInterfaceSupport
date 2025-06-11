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
  requires(data_size > 0) constexpr bool FunctionDataSize<data_size>::isScalar()
      const noexcept {
    return data_size == 1;
  }  // end of isScalar

  template <size_type data_size>
  requires(data_size > 0) constexpr size_type
      FunctionDataSize<data_size>::getNumberOfComponents() const noexcept {
    return data_size;
  }  // end of FunctionDataSize<dynamic_extent>::getNumberOfComponents

  constexpr bool FunctionDataSize<dynamic_extent>::isScalar() const noexcept {
    return this->getNumberOfComponents() == 1;
  }  // end of isScalar

  constexpr size_type FunctionDataSize<dynamic_extent>::getNumberOfComponents()
      const noexcept {
    return this->data_size;
  }  // end of FunctionDataSize<dynamic_extent>::getNumberOfComponents

  template <size_type data_stride>
  requires(data_stride > 0) constexpr size_type
      FunctionDataStride<data_stride>::getDataStride() const noexcept {
    return data_stride;
  }  // end of getDataStride

  constexpr size_type FunctionDataStride<dynamic_extent>::getDataStride()
      const noexcept {
    return this->data_stride;
  }  // end of getDataStride

  constexpr bool has_dynamic_properties(const DataLayoutDescription& layout) {
    return (layout.data_size == dynamic_extent) ||  //
           (layout.data_stride == dynamic_extent);
  }  // end of has_dynamic_properties

  template <DataLayoutDescription layout>
  requires((layout.data_size > 0) && (layout.data_stride > 0))  //
      constexpr size_type
      DataLayout<layout>::getDataOffset(const size_type i) const noexcept {
    if constexpr (layout.data_stride != dynamic_extent) {
      return i * layout.data_stride;
    } else {
      return i * (this->getDataStride());
    }
  }  // end of getDataOffset

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr bool FunctionView<Space, layout, is_mutable>::
          checkPreconditions(AbstractErrorHandler& eh,
                             const Space& s,
                             typename FunctionView::ExternalData v)  //
      requires((layout.data_size != dynamic_extent) &&
               (layout.data_stride != dynamic_extent)) {
    if (s.get() == nullptr) {
      return eh.registerErrorMessage("invalid space");
    }
    const auto space_size = getSpaceSize(*s);
    if (space_size == 0) {
      // this may happen due to partionning in parallel
      return true;
    }
    return v.size() >= layout.data_stride * (space_size - 1) + layout.data_size;
  }

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr FunctionView<Space, layout, is_mutable>::FunctionView(
          const Space& s,
          typename FunctionView::ExternalData v)  //
      requires((layout.data_size != dynamic_extent) &&
               (layout.data_stride != dynamic_extent))
      : FunctionView(preconditions_check, s, v) {}

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      template <bool doPreconditionsCheck>
      constexpr FunctionView<Space, layout, is_mutable>::FunctionView(
          const PreconditionsCheck<doPreconditionsCheck>& pcheck,
          const Space& s,
          typename FunctionView::ExternalData v)  //
      requires((layout.data_size != dynamic_extent) &&
               (layout.data_stride != dynamic_extent))
      : PreconditionsChecker<FunctionView>(pcheck, s, v), space(s), values(v) {}

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr bool FunctionView<Space, layout, is_mutable>::
          checkPreconditions(AbstractErrorHandler& eh,
                             const Space& s,
                             FunctionView::ExternalData v,
                             const size_type dsize)  //
      requires((layout.data_size == dynamic_extent) &&
               (layout.data_stride != dynamic_extent)) {
    if (s.get() == nullptr) {
      return eh.registerErrorMessage("invalid space");
    }
    if (dsize <= 0) {
      return eh.registerErrorMessage("invalid number of components");
    }
    const auto space_size = s->size();
    if (space_size == 0) {
      // this may happen due to partionning in parallel
      return true;
    }
    const auto min_size = (layout.data_stride) * (space_size - 1) + dsize;
    return v.size() >= min_size;
  }  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr FunctionView<Space, layout, is_mutable>::FunctionView(
          const Space& s,
          FunctionView::ExternalData v,
          const size_type dsize)  //
      requires((layout.data_size == dynamic_extent) &&
               (layout.data_stride != dynamic_extent))
      : FunctionView(preconditions_check, s, v, dsize) {}

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      template <bool doPreconditionsCheck>
      constexpr FunctionView<Space, layout, is_mutable>::FunctionView(
          const PreconditionsCheck<doPreconditionsCheck>& pcheck,
          const Space& s,
          FunctionView::ExternalData v,
          const size_type dsize)  //
      requires((layout.data_size == dynamic_extent) &&
               (layout.data_stride != dynamic_extent))
      : PreconditionsChecker<FunctionView>(pcheck, s, dsize),
        space(s),
        values(v) {
    this->data_size = dsize;
  }  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr bool FunctionView<Space, layout, is_mutable>::
          checkPreconditions(AbstractErrorHandler& eh,
                             const Space& s,
                             FunctionView::ExternalData v,
                             const size_type dstride)  //
      requires((layout.data_size != dynamic_extent) &&
               (layout.data_stride == dynamic_extent)) {
    if (s.get() == nullptr) {
      return eh.registerErrorMessage("invalid space");
    }
    if (dstride <= 0) {
      return eh.registerErrorMessage("invalid stride");
    }
    const auto space_size = getSpaceSize(s);
    if (space_size == 0) {
      // this may happen due to partionning in parallel
      return true;
    }
    const auto min_size = dstride * (space_size - 1) + layout.data_size;
    if (v.size() < min_size) {
      return eh.registerErrorMessage("invalid external data size");
    }
    return true;
  }  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr FunctionView<Space, layout, is_mutable>::FunctionView(
          const Space& s,
          FunctionView::ExternalData v,
          const size_type dstride)  //
      requires((layout.data_size != dynamic_extent) &&
               (layout.data_stride == dynamic_extent))
      : FunctionView(preconditions_check, s, v, dstride) {}

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      template <bool doPreconditionsCheck>
      constexpr FunctionView<Space, layout, is_mutable>::FunctionView(
          const PreconditionsCheck<doPreconditionsCheck>& pcheck,
          const Space& s,
          FunctionView::ExternalData v,
          const size_type dstride)  //
      requires((layout.data_size != dynamic_extent) &&
               (layout.data_stride == dynamic_extent))
      : PreconditionsChecker<FunctionView>(pcheck, s, v, dstride),
        space(s),
        values(v) {
    this->data_stride = dstride;
  }  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr bool FunctionView<Space, layout, is_mutable>::
          checkPreconditions(AbstractErrorHandler& eh,
                             const Space& s,
                             FunctionView::ExternalData v,
                             const size_type dsize,
                             const size_type dstride)  //
      requires((layout.data_size == dynamic_extent) &&
               (layout.data_stride == dynamic_extent)) {
    if (dsize <= 0) {
      return eh.registerErrorMessage("invalid number of components");
    }
    if (dstride <= 0) {
      return eh.registerErrorMessage("invalid stride");
    }
    if (dsize > dstride) {
      return eh.registerErrorMessage(
          "the number of components is greater than the stride");
    }
    const auto space_size = getSpaceSize(s);
    if (space_size == 0) {
      // this may happen due to partionning in parallel
      return true;
    }
    const auto min_size = dstride * (space_size - 1) + dsize;
    if (v.size() < min_size) {
      return eh.registerErrorMessage("invalid external data size");
    }
    return true;
  }  // end of checkPreconditions

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr bool FunctionView<Space, layout, is_mutable>::
          checkPreconditions(AbstractErrorHandler& eh,
                             const Space& s,
                             FunctionView::ExternalData v,
                             const size_type dsize)  //
      requires((layout.data_size == dynamic_extent) &&
               (layout.data_stride == dynamic_extent)) {
    return checkPreconditions(eh, s, v, dsize, dsize);
  }  // end of checkPreconditions

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr FunctionView<Space, layout, is_mutable>::FunctionView(
          const Space& s,
          FunctionView::ExternalData v,
          const size_type dsize,
          const size_type dstride)  //
      requires((layout.data_size == dynamic_extent) &&
               (layout.data_stride == dynamic_extent))
      : FunctionView(preconditions_check, s, v, dsize, dstride) {}

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      template <bool doPreconditionsCheck>
      constexpr FunctionView<Space, layout, is_mutable>::FunctionView(
          const PreconditionsCheck<doPreconditionsCheck>& pcheck,
          const Space& s,
          FunctionView::ExternalData v,
          const size_type dsize,
          const size_type dstride)  //
      requires((layout.data_size == dynamic_extent) &&
               (layout.data_stride == dynamic_extent))
      : PreconditionsChecker<FunctionView>(pcheck, s, v, dstride, dstride),
        space(s),
        values(v) {
    this->data_size = dsize;
    this->data_stride = dstride;
  }  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr FunctionView<Space, layout, is_mutable>::FunctionView(
          const Space& s,
          FunctionView::ExternalData v,
          const size_type dsize)  //
      requires((layout.data_size == dynamic_extent) &&
               (layout.data_stride == dynamic_extent))
      : FunctionView(preconditions_check, s, v, dsize, dsize) {
  }  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      template <bool doPreconditionsCheck>
      constexpr FunctionView<Space, layout, is_mutable>::FunctionView(
          const PreconditionsCheck<doPreconditionsCheck>& pcheck,
          const Space& s,
          FunctionView::ExternalData v,
          const size_type dsize)  //
      requires((layout.data_size == dynamic_extent) &&
               (layout.data_stride == dynamic_extent))
      : FunctionView(pcheck, s, v, dsize, dsize) {}  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr FunctionView<Space, layout, is_mutable>::FunctionView(
          const Space& s,
          FunctionView::ExternalData v,
          const DataLayout<layout>& l)  //
      requires((layout.data_size == dynamic_extent) &&
               (layout.data_stride == dynamic_extent))
      : FunctionView(preconditions_check, s, v, l.size, l.stride) {
  }  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      template <bool doPreconditionsCheck>
      constexpr FunctionView<Space, layout, is_mutable>::FunctionView(
          const PreconditionsCheck<doPreconditionsCheck>& pcheck,
          const Space& s,
          FunctionView::ExternalData v,
          const DataLayout<layout>& l)  //
      requires((layout.data_size == dynamic_extent) &&
               (layout.data_stride == dynamic_extent))
      : FunctionView(pcheck, s, v, l.size, l.stride) {}  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr FunctionView<Space, layout, is_mutable>::FunctionView(
          const Space& s,
          FunctionView::ExternalData v,
          const DataLayout<layout>& l)  //
      requires((layout.data_size != dynamic_extent) &&
               (layout.data_stride == dynamic_extent))
      : FunctionView(preconditions_check, s, v, l.stride) {}

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      template <bool doPreconditionsCheck>
      constexpr FunctionView<Space, layout, is_mutable>::FunctionView(
          const PreconditionsCheck<doPreconditionsCheck>& pcheck,
          const Space& s,
          FunctionView::ExternalData v,
          const DataLayout<layout>& l)  //
      requires((layout.data_size != dynamic_extent) &&
               (layout.data_stride == dynamic_extent))
      : FunctionView(pcheck, s, v, l.stride) {}

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr FunctionView<Space, layout, is_mutable>::FunctionView(
          const Space& s,
          FunctionView::ExternalData v,
          const DataLayout<layout>& l)  //
      requires((layout.data_size == dynamic_extent) &&
               (layout.data_stride != dynamic_extent))
      : FunctionView(preconditions_check, s, v, l.size) {
  }  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      template <bool doPreconditionsCheck>
      constexpr FunctionView<Space, layout, is_mutable>::FunctionView(
          const PreconditionsCheck<doPreconditionsCheck>& pcheck,
          const Space& s,
          FunctionView::ExternalData v,
          const DataLayout<layout>& l)  //
      requires((layout.data_size == dynamic_extent) &&
               (layout.data_stride != dynamic_extent))
      : FunctionView(pcheck, s, v, l.size) {}  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr bool FunctionView<Space, layout, is_mutable>::
          checkCompatibility(const FunctionView& v) const {
    if (areEquivalent(this->getSpace(), v.getSpacePointer())) {
      return false;
    }
    return this->data_size == v.getNumberOfComponents();
  }  // end of checkCompatibility

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr const Space& FunctionView<Space, layout, is_mutable>::getSpace()
          const noexcept {
    return this->space;
  }  // end of getSpace

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      bool FunctionView<Space, layout, is_mutable>::check(
          Context&) const noexcept {
    return true;
  }  // end of allocateWorkspace

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr void FunctionView<Space, layout, is_mutable>::
          allocateWorkspace() noexcept {}  // end of allocateWorkspace

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr real* FunctionView<Space, layout, is_mutable>::data(
          mgis::attributes::UnsafeAttribute,
          const size_type o)  //
      requires(is_mutable&& LinearElementSpaceConcept<Space> &&
               (!hasElementWorkspace<Space>)) {
    return this->values.data() + this->getDataOffset(o);
  }  // end of data

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr real* FunctionView<Space, layout, is_mutable>::data(
          mgis::attributes::UnsafeAttribute,
          const size_type e,
          const size_type i)  //
      requires(is_mutable&& LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return this->values.data() +
           this->getDataOffset(this->space->getQuadraturePointOffset(e, i));
  }  // end of data

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr const real* FunctionView<Space, layout, is_mutable>::data(
          mgis::attributes::UnsafeAttribute,
          const size_type o) const  //
      requires(LinearElementSpaceConcept<Space> &&
               (!hasElementWorkspace<Space>)) {
    return this->values.data() + this->getDataOffset(o);
  }  // end of data

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr const real* FunctionView<Space, layout, is_mutable>::data(
          mgis::attributes::UnsafeAttribute,
          const size_type e,
          const size_type i) const
      requires(LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return this->values.data() +
           this->getDataOffset(getQuadraturePointOffset(this->space, e, i));
  }  // end of data

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr typename FunctionView<Space, layout, is_mutable>::ValuesView
      FunctionView<Space, layout, is_mutable>::operator()(const size_type o)  //
      requires(is_mutable&& LinearElementSpaceConcept<Space> &&
               (!hasElementWorkspace<Space>)) {
    if constexpr (layout.data_size == 1) {
      return *(this->values.data() + this->getDataOffset(o));
    } else {
      return ValuesView(this->values.data() + this->getDataOffset(o),
                        this->getNumberOfComponents());
    }
  }  // end of operator()

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      typename FunctionView<Space, layout, is_mutable>::ValuesView
      constexpr FunctionView<Space, layout, is_mutable>::operator()(
          const size_type e,
          const size_type i)  //
      requires(is_mutable&& LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return this->operator()(this->space->getQuadraturePointOffset(e, i));
  }  // end of operator()

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      typename FunctionView<Space, layout, is_mutable>::ConstValuesView
      constexpr FunctionView<Space, layout, is_mutable>::operator()(
          const size_type o) const requires(LinearElementSpaceConcept<Space> &&
                                            (!hasElementWorkspace<Space>)) {
    if constexpr (layout.data_size == 1) {
      return *(this->values.data() + this->getDataOffset(o));
    } else {
      return ConstValuesView(this->values.data() + this->getDataOffset(o),
                             this->getNumberOfComponents());
    }
  }  // end of operator()

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      typename FunctionView<Space, layout, is_mutable>::ConstValuesView
      constexpr FunctionView<Space, layout, is_mutable>::operator()(
          const size_type e, const size_type i) const
      requires(LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return this->operator()(this->space->getQuadraturePointOffset(e, i));
  }  // end of operator()

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      constexpr std::span<const real>            //
      FunctionView<Space, layout, is_mutable>::data() const {
    return this->values;
  }

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      constexpr Function<Space, N>::Function(const Space& s) requires(
          N != dynamic_extent)
      : Function(preconditions_check, s) {}  // end of Function

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      constexpr Function<Space, N>::Function(
          const Space& s, const size_type dsize) requires(N == dynamic_extent)
      : Function(preconditions_check, s, dsize) {}  // end of Function

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      constexpr bool Function<Space, N>::checkPreconditions(
          AbstractErrorHandler&, const Space&) requires(N != dynamic_extent) {
    return true;
  }

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      template <bool doPreconditionsCheck>
      constexpr Function<Space, N>::Function(
          const PreconditionsCheck<doPreconditionsCheck>& pcheck,
          const Space& s)  //
      requires(N != dynamic_extent)
      : PreconditionsChecker<Function>(pcheck, s),
        FunctionStorage<Space, N>(s),
        FunctionView<Space, {.data_size = N, .data_stride = N}, true>(
            pcheck, s, this->storage_values) {}

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      constexpr bool Function<Space, N>::checkPreconditions(
          AbstractErrorHandler& eh,
          const Space&,
          const size_type dsize) requires(N == dynamic_extent) {
    if (dsize <= 0) {
      return eh.registerErrorMessage("invalid number of components");
    }
    return true;
  }

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      template <bool doPreconditionsCheck>
      constexpr Function<Space, N>::Function(
          const PreconditionsCheck<doPreconditionsCheck>& pcheck,
          const Space& s,
          const size_type dsize)  //
      requires(N == dynamic_extent)
      : PreconditionsChecker<Function>(pcheck, s, dsize),
        FunctionStorage<Space, N>(s, dsize),
        FunctionView<Space, {.data_size = N, .data_stride = N}, true>(
            pcheck, s, this->storage_values, dsize, dsize) {}

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      constexpr Function<Space, N>::Function(const Function& f) requires(
          N == dynamic_extent)
      : PreconditionsChecker<Function>(f),
        FunctionStorage<Space, N>(f),
        FunctionView<Space, {.data_size = N, .data_stride = N}, true>(
            f.getSpace(),
            this->storage_values,
            f.getNumberOfComponents(),
            f.getNumberOfComponents()) {}  // end of Function

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      constexpr Function<Space, N>::Function(Function&& f) requires(
          N == dynamic_extent)
      : PreconditionsChecker<Function>(f),
        FunctionStorage<Space, N>(std::forward<FunctionStorage<Space, N>>(f)),
        FunctionView<Space, {.data_size = N, .data_stride = N}, true>(
            f.getSpace(),
            this->storage_values,
            f.getNumberOfComponents(),
            f.getNumberOfComponents()) {}  // end of Function

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      constexpr Function<Space, N>::Function(const Function& f) requires(
          N != dynamic_extent)
      : PreconditionsChecker<Function>(f),
        FunctionStorage<Space, N>(f),
        FunctionView<Space, {.data_size = N, .data_stride = N}, true>(
            no_precondition_check, f.getSpace(), this->storage_values) {
  }  // end of Function

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      constexpr Function<Space, N>::Function(Function&& f) requires(
          N != dynamic_extent)
      : PreconditionsChecker<Function>(f),
        FunctionStorage<Space, N>(std::forward<FunctionStorage<Space, N>>(f)),
        FunctionView<Space, {.data_size = N, .data_stride = N}, true>(
            no_precondition_check, f.getSpace(), this->storage_values) {
  }  // end of Function

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      constexpr FunctionView<Space,
                             {.data_size = N, .data_stride = N},
                             true>  //
      Function<Space, N>::view() {
    return *this;
  }  // end of view

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      constexpr FunctionView<Space,
                             {.data_size = N, .data_stride = N},
                             false>  //
      Function<Space, N>::view() const {
    if constexpr (N == dynamic_extent) {
      return FunctionView<Space, {.data_size = N, .data_stride = N}, false>(
          this->getSpace(), this->storage_values, this->getNumberOfComponents(),
          this->getNumberOfComponents());
    } else {
      return FunctionView<Space, {.data_size = N, .data_stride = N}, false>(
          this->getSpace(), this->storage_values);
    }
  }  // end of view

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      constexpr Function<Space, N>::~Function() = default;

  template <FunctionalSpaceConcept Space, size_type N>
  constexpr auto view(const Function<Space, N>& f) {
    return f.view();
  }  // end of view

  template <size_type N, FunctionalSpaceConcept Space, size_type N2>
  constexpr auto view(const Function<Space, N2>& f)  //
      requires((N > 0) && (N != dynamic_extent) && (N == N2)) {
    return f.view();
  }  // end of view

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  constexpr const auto& getSpace(
      const FunctionView<Space, layout, is_mutable>& f) {
    return f.getSpace();
  }  // end of getSpace

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_FUNCTION_IXX */
