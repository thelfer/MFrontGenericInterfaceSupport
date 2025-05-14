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

  inline bool FunctionDataSize<dynamic_extent>::isScalar() const noexcept {
    return this->getNumberOfComponents() == 1;
  }  // end of isScalar

  inline size_type FunctionDataSize<dynamic_extent>::getNumberOfComponents()
      const noexcept {
    return this->data_size;
  }  // end of FunctionDataSize<dynamic_extent>::getNumberOfComponents

  template <size_type data_stride>
  requires(data_stride > 0) constexpr size_type
      FunctionDataStride<data_stride>::getDataStride() const noexcept {
    return data_stride;
  }  // end of getDataStride

  inline size_type FunctionDataStride<dynamic_extent>::getDataStride()
      const noexcept {
    return this->data_stride;
  }  // end of getDataStride

  constexpr bool has_dynamic_properties(const DataLayoutDescription& layout) {
    return (layout.data_size == dynamic_extent) ||  //
           (layout.data_stride == dynamic_extent);
  }  // end of has_dynamic_properties

  template <DataLayoutDescription layout>
  requires((layout.data_size > 0) && (layout.data_stride > 0)) inline size_type
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
      bool FunctionView<Space, layout, is_mutable>::checkPreconditions(
          Context& ctx,
          std::shared_ptr<const Space> s,
          typename FunctionView::ExternalData v) requires((layout.data_size !=
                                                           dynamic_extent) &&
                                                          (layout.data_stride !=
                                                           dynamic_extent)) {
    if (s.get() == nullptr) {
      return ctx.registerErrorMessage("invalid space");
    }
    const auto space_size = s->size();
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
      FunctionView<Space, layout, is_mutable>::FunctionView(
          std::shared_ptr<const Space> s, typename FunctionView::ExternalData v)
  requires((layout.data_size != dynamic_extent) &&
           (layout.data_stride != dynamic_extent))
      : space(s), values(v) {
    Context ctx;
    if (!checkPreconditions(ctx, s, v)) {
      raise(ctx.getErrorMessage());
    }
  }

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      bool FunctionView<Space, layout, is_mutable>::checkPreconditions(
          Context& ctx,
          std::shared_ptr<const Space> s,
          FunctionView::ExternalData v,
          const size_type dsize) requires((layout.data_size ==
                                           dynamic_extent) &&
                                          (layout.data_stride !=
                                           dynamic_extent)) {
    if (s.get() == nullptr) {
      return ctx.registerErrorMessage("invalid space");
    }
    if (dsize <= 0) {
      return ctx.registerErrorMessage("invalid number of components");
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
      FunctionView<Space, layout, is_mutable>::FunctionView(
          std::shared_ptr<const Space> s,
          FunctionView::ExternalData v,
          const size_type dsize)
  requires((layout.data_size == dynamic_extent) &&
           (layout.data_stride != dynamic_extent))
      : space(s), values(v) {
    Context ctx;
    if (!checkPreconditions(ctx, s, v, dsize)) {
      raise(ctx.getErrorMessage());
    }
    this->data_size = dsize;
  }  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      bool FunctionView<Space, layout, is_mutable>::checkPreconditions(
          Context& ctx,
          std::shared_ptr<const Space> s,
          FunctionView::ExternalData v,
          const size_type dstride) requires((layout.data_size !=
                                             dynamic_extent) &&
                                            (layout.data_stride ==
                                             dynamic_extent)) {
    if (s.get() == nullptr) {
      return ctx.registerErrorMessage("invalid space");
    }
    if (dstride <= 0) {
      return ctx.registerErrorMessage("invalid stride");
    }
    const auto space_size = s->size();
    if (space_size == 0) {
      // this may happen due to partionning in parallel
      return true;
    }
    const auto min_size = dstride * (space_size - 1) + layout.data_size;
    if (v.size() < min_size) {
      return ctx.registerErrorMessage("invalid external data size");
    }
    return true;
  }  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      FunctionView<Space, layout, is_mutable>::FunctionView(
          std::shared_ptr<const Space> s,
          FunctionView::ExternalData v,
          const size_type dstride)
  requires((layout.data_size != dynamic_extent) &&
           (layout.data_stride == dynamic_extent))
      : space(s), values(v) {
    Context ctx;
    if (!checkPreconditions(ctx, s, v, dstride)) {
      raise(ctx.getErrorMessage());
    }
    this->data_stride = dstride;
  }  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      bool FunctionView<Space, layout, is_mutable>::checkPreconditions(
          Context& ctx,
          std::shared_ptr<const Space> s,
          FunctionView::ExternalData v,
          const size_type dsize,
          const size_type dstride) requires((layout.data_size ==
                                             dynamic_extent) &&
                                            (layout.data_stride ==
                                             dynamic_extent)) {
    if (s.get() == nullptr) {
      return ctx.registerErrorMessage("invalid space");
    }
    if (dsize <= 0) {
      return ctx.registerErrorMessage("invalid number of components");
    }
    if (dstride <= 0) {
      return ctx.registerErrorMessage("invalid stride");
    }
    if (dsize > dstride) {
      return ctx.registerErrorMessage(
          "the number of components is greater than the stride");
    }
    const auto space_size = s->size();
    if (space_size == 0) {
      // this may happen due to partionning in parallel
      return true;
    }
    const auto min_size = dstride * (space_size - 1) + dsize;
    if (v.size() < min_size) {
      return ctx.registerErrorMessage("invalid external data size");
    }
    return true;
  }  // end of checkPreconditions

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      bool FunctionView<Space, layout, is_mutable>::checkPreconditions(
          Context& ctx,
          std::shared_ptr<const Space> s,
          FunctionView::ExternalData v,
          const size_type dsize) requires((layout.data_size ==
                                           dynamic_extent) &&
                                          (layout.data_stride ==
                                           dynamic_extent)) {
    return checkPreconditions(ctx, s, v, dsize, dsize);
  }  // end of checkPreconditions

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      FunctionView<Space, layout, is_mutable>::FunctionView(
          std::shared_ptr<const Space> s,
          FunctionView::ExternalData v,
          const size_type dsize,
          const size_type dstride)
  requires((layout.data_size == dynamic_extent) &&
           (layout.data_stride == dynamic_extent))
      : space(s), values(v) {
    Context ctx;
    if (!checkPreconditions(ctx, s, v, dsize, dstride)) {
      raise(ctx.getErrorMessage());
    }
    this->data_size = dsize;
    this->data_stride = dstride;
  }  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      FunctionView<Space, layout, is_mutable>::FunctionView(
          std::shared_ptr<const Space> s,
          FunctionView::ExternalData v,
          const size_type dsize)
  requires((layout.data_size == dynamic_extent) &&
           (layout.data_stride == dynamic_extent))
      : FunctionView(s, v, dsize, dsize) {}  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      FunctionView<Space, layout, is_mutable>::FunctionView(
          std::shared_ptr<const Space> s,
          FunctionView::ExternalData v,
          const DataLayout<layout>& l)
  requires((layout.data_size == dynamic_extent) &&
           (layout.data_stride == dynamic_extent))
      : FunctionView(s, v, l.size, l.stride) {}  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      FunctionView<Space, layout, is_mutable>::FunctionView(
          std::shared_ptr<const Space> s,
          FunctionView::ExternalData v,
          const DataLayout<layout>& l)
  requires((layout.data_size != dynamic_extent) &&
           (layout.data_stride == dynamic_extent))
      : FunctionView(s, v, l.stride) {}

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      FunctionView<Space, layout, is_mutable>::FunctionView(
          std::shared_ptr<const Space> s,
          FunctionView::ExternalData v,
          const DataLayout<layout>& l)
  requires((layout.data_size == dynamic_extent) &&
           (layout.data_stride != dynamic_extent))
      : FunctionView(s, v, l.size) {}  // end of FunctionView

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
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
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      const Space& FunctionView<Space, layout, is_mutable>::getSpace()
          const noexcept {
    return *(this->space);
  }  // end of getSpace

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      std::shared_ptr<const Space> FunctionView<Space, layout, is_mutable>::
          getSpacePointer()
  const noexcept { return this->space; }  // end of getSpacePointer

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
      void FunctionView<Space, layout, is_mutable>::
          allocateWorkspace() noexcept {}  // end of allocateWorkspace

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      real& FunctionView<Space, layout, is_mutable>::getValue(
          const size_type o) requires(allowScalarAccessor&& is_mutable&&
                                          LinearElementSpaceConcept<Space> &&
                                      (!hasElementWorkspace<Space>)) {
    return *(this->values.data() + this->getDataOffset(o));
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      typename FunctionView<Space, layout, is_mutable>::ValuesView
      FunctionView<Space, layout, is_mutable>::getValues(
          const size_type
              o) requires(is_mutable&& LinearElementSpaceConcept<Space> &&
                          (!hasElementWorkspace<Space>)) {
    if constexpr (layout.data_size == 1) {
      return *(this->values.data() + this->getDataOffset(o));
    } else {
      return ValuesView(this->values.data() + this->getDataOffset(o),
                        this->getNumberOfComponents());
    }
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      template <size_type N>
      std::span<real, N> FunctionView<Space, layout, is_mutable>::getValues(
          const size_type o) requires((layout.data_size == dynamic_extent) &&
                                      is_mutable &&
                                      LinearElementSpaceConcept<Space> &&
                                      (!hasElementWorkspace<Space>)) {
    return std::span<real, N>(this->values.data() + this->getDataOffset(o),
                              this->getNumberOfComponents());
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
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
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
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
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      template <size_type N>
      std::span<real, N> FunctionView<Space, layout, is_mutable>::getValues(
          const size_type e,
          const size_type i) requires((layout.data_size == dynamic_extent) &&
                                      is_mutable &&
                                      LinearQuadratureSpaceConcept<Space> &&
                                      (!hasCellWorkspace<Space>)) {
    return this->template getValues<N>(
        this->space->getQuadraturePointOffset(e, i));
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      const real& FunctionView<Space, layout, is_mutable>::getValue(
          const size_type o) const
      requires(allowScalarAccessor&& LinearElementSpaceConcept<Space> &&
               (!hasElementWorkspace<Space>)) {
    return *(this->values.data() + this->getDataOffset(o));
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      typename FunctionView<Space, layout, is_mutable>::ConstValuesView
      FunctionView<Space, layout, is_mutable>::getValues(
          const size_type o) const requires(LinearElementSpaceConcept<Space> &&
                                            (!hasElementWorkspace<Space>)) {
    if constexpr (layout.data_size == 1) {
      return *(this->values.data() + this->getDataOffset(o));
    } else {
      return ConstValuesView(this->values.data() + this->getDataOffset(o),
                             this->getNumberOfComponents());
    }
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      template <size_type N>
      std::span<const real, N> FunctionView<Space, layout, is_mutable>::
          getValues(const size_type o) const
      requires(LinearElementSpaceConcept<Space> &&
               (!hasElementWorkspace<Space>)) {
    return std::span<const real, N>(
        this->values.data() + this->getDataOffset(o),
        this->getNumberOfComponents());
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      const real& FunctionView<Space, layout, is_mutable>::getValue(
          const size_type e, const size_type i) const
      requires(allowScalarAccessor&& LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return this->getValue(this->space->getQuadraturePointOffset(e, i));
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      typename FunctionView<Space, layout, is_mutable>::ConstValuesView
      FunctionView<Space, layout, is_mutable>::getValues(
          const size_type e, const size_type i) const
      requires(LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return this->getValues(this->space->getQuadraturePointOffset(e, i));
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      template <size_type N>
      std::span<const real, N> FunctionView<Space, layout, is_mutable>::
          getValues(const size_type e, const size_type i) const
      requires(LinearQuadratureSpaceConcept<Space> &&
               (!hasCellWorkspace<Space>)) {
    return this->template getValues<N>(
        this->space->getQuadraturePointOffset(e, i));
  }  // end of getValues

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      typename FunctionView<Space, layout, is_mutable>::ValuesView
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
      FunctionView<Space, layout, is_mutable>::operator()(const size_type e,
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
      FunctionView<Space, layout, is_mutable>::operator()(
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
      FunctionView<Space, layout, is_mutable>::operator()(
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
      std::span<const real> FunctionView<Space, layout, is_mutable>::getValues()
  const { return this->values; }

  template <FunctionalSpaceConcept Space,
            DataLayoutDescription layout,
            bool is_mutable>
  requires(LinearElementSpaceConcept<Space> ||
           LinearQuadratureSpaceConcept<Space>)  //
      FunctionView<Space, layout, is_mutable>::~FunctionView() noexcept =
          default;

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      bool Function<Space, N>::checkPreconditions(
          Context& ctx,
          const std::shared_ptr<const Space>& s) requires(N != dynamic_extent) {
    if (isInvalid(s)) {
      return ctx.registerErrorMessage("invalid space pointer");
    }
    return true;
  }

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))        //
      Function<Space, N>::Function(std::shared_ptr<const Space> s)  //
  requires(N != dynamic_extent)
      : FunctionStorage<Space, N>(*s),
        FunctionView<Space, {.data_size = N, .data_stride = N}, true>(
            s, this->storage_values) {
    Context ctx;
    if (!checkPreconditions(ctx, s)) {
      raise(ctx.getErrorMessage());
    }
  }

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      bool Function<Space, N>::checkPreconditions(
          Context& ctx,
          const std::shared_ptr<const Space>& s,
          const size_type dsize) requires(N == dynamic_extent) {
    if (isInvalid(s)) {
      return ctx.registerErrorMessage("invalid space pointer");
    }
    if (dsize <= 0) {
      return ctx.registerErrorMessage("invalid number of components");
    }
    return true;
  }

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      Function<Space, N>::Function(std::shared_ptr<const Space> s,
                                   const size_type dsize)  //
  requires(N == dynamic_extent)
      : FunctionStorage<Space, N>(*s, dsize),
        FunctionView<Space, {.data_size = N, .data_stride = N}, true>(
            s, this->storage_values, dsize, dsize) {
    Context ctx;
    if (!checkPreconditions(ctx, s, dsize)) {
      raise(ctx.getErrorMessage());
    }
  }

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      Function<Space, N>::Function(const Function& f)
  requires(N == dynamic_extent)
      : FunctionStorage<Space, N>(f),
        FunctionView<Space, {.data_size = N, .data_stride = N}, true>(
            f.getSpacePointer(),
            this->storage_values,
            f.getNumberOfComponents(),
            f.getNumberOfComponents()) {}  // end of Function

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      Function<Space, N>::Function(Function&& f)
  requires(N == dynamic_extent)
      : FunctionStorage<Space, N>(std::forward<FunctionStorage<Space, N>>(f)),
        FunctionView<Space, {.data_size = N, .data_stride = N}, true>(
            f.getSpacePointer(),
            this->storage_values,
            f.getNumberOfComponents(),
            f.getNumberOfComponents()) {}  // end of Function

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      Function<Space, N>::Function(const Function& f)
  requires(N != dynamic_extent)
      : FunctionStorage<Space, N>(f),
        FunctionView<Space, {.data_size = N, .data_stride = N}, true>(
            f.getSpacePointer(), this->storage_values) {}  // end of Function

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      Function<Space, N>::Function(Function&& f)
  requires(N != dynamic_extent)
      : FunctionStorage<Space, N>(std::forward<FunctionStorage<Space, N>>(f)),
        FunctionView<Space, {.data_size = N, .data_stride = N}, true>(
            f.getSpacePointer(), this->storage_values) {}  // end of Function

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      FunctionView<Space,
                   {.data_size = N, .data_stride = N},
                   true> Function<Space, N>::view() {
    return *this;
  }  // end of view

  template <FunctionalSpaceConcept Space, size_type N>
  requires((N > 0) && (LinearElementSpaceConcept<Space> ||
                       LinearQuadratureSpaceConcept<Space>))  //
      FunctionView<Space,
                   {.data_size = N, .data_stride = N},
                   false> Function<Space, N>::view()
  const {
    if constexpr (N == dynamic_extent) {
      return FunctionView<Space, {.data_size = N, .data_stride = N}, false>(
          this->getSpacePointer(), this->storage_values,
          this->getNumberOfComponents(), this->getNumberOfComponents());
    } else {
      return FunctionView<Space, {.data_size = N, .data_stride = N}, false>(
          this->getSpacePointer(), this->storage_values);
    }
  }  // end of view

  template <FunctionalSpaceConcept Space, size_type N>
  auto view(const Function<Space, N>& f) {
    return f.view();
  }  // end of view

  template <size_type N, FunctionalSpaceConcept Space, size_type N2>
  auto view(const Function<Space, N2>& f)  //
      requires((N > 0) && (N != dynamic_extent) && (N == N2)) {
    return f.view();
  }  // end of view

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_FUNCTION_IXX */
