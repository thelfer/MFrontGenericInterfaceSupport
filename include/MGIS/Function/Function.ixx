/*!
 * \file   MGIS/Function/Function.ixx
 * \brief
 * \author Thomas Helfer
 * \date   10/03/2021
 */

#ifndef LIB_MGIS_FUNCTION_FUNCTION_IXX
#define LIB_MGIS_FUNCTION_FUNCTION_IXX

namespace mgis::function {

  template <size_type offset>
  constexpr size_type FunctionDataOffset<offset>::getDataOffset()
      const noexcept {
    return offset;
  }  // end of getDataOffset

  inline size_type FunctionDataOffset<dynamic_extent>::getDataOffset()
      const noexcept {
    return this->data_begin;
  }  // end of getDataOffset

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
    return this->data_stride;
  }  // end of getDataStride

  inline size_type FunctionDataStride<dynamic_extent>::getDataStride()
      const noexcept {
    return this->data_stride;
  }  // end of getDataStride

  constexpr bool has_dynamic_properties(const DataLayoutDescription& layout) {
    return (layout.size == dynamic_extent) ||    //
           (layout.offset == dynamic_extent) ||  //
           (layout.stride == dynamic_extent);
  }  // end of has_dynamic_properties

  template <DataLayoutDescription layout>
  inline size_type DataLayout<layout>::getDataOffset(
      const size_type i) const noexcept {
    if constexpr ((layout.stride != dynamic_extent) &&
                  (layout.offset != dynamic_extent)) {
      if constexpr (layout.offset == 0) {
        return i * layout.stride;
      } else {
        return i * layout.stride + layout.offset;
      }
    } else if constexpr (layout.offset != dynamic_extent) {
      if constexpr (layout.offset == 0) {
        return i * (this->getDataStride());
      } else {
        return i * (this->getDataStride()) + layout.offset;
      }
    } else {
      return i * (this->getDataStride()) + this->getDataOffset();
    }
  }  // end of getDataOffset

  inline const AbstractSpace& ImmutableFunctionView::getSpace() const noexcept {
    return *(this->qspace);
  }  // end of getSpace

  inline std::shared_ptr<const AbstractSpace>
  ImmutableFunctionView::getSpacePointer() const noexcept {
    return this->qspace;
  }  // end of getSpacePointer

  inline const real& ImmutableFunctionView::getValue(const size_type o) const {
    return *(this->immutable_values.data() + this->getDataOffset(o));
  }  // end of getValues

  inline std::span<const real> ImmutableFunctionView::getValues(
      const size_type o) const {
    return std::span<const real>(
        this->immutable_values.data() + this->getDataOffset(o),
        this->getNumberOfComponents());
  }  // end of getValues

  template <size_type N>
  inline std::span<const real, N> ImmutableFunctionView::getValues(
      const size_type o) const {
    return std::span<const real, N>(
        this->immutable_values.data() + this->getDataOffset(o),
        this->getNumberOfComponents());
  }  // end of getValues

  inline std::span<const real> ImmutableFunctionView::getValues() const {
    return this->immutable_values;
  }

  inline real& Function::getValue(const size_type o) {
    return *(this->values.data() + this->getDataOffset(o));
  }  // end of getValues

  inline std::span<real> Function::getValues(const size_type o) {
    return std::span<real>(this->values.data() + this->getDataOffset(o),
                           this->getNumberOfComponents());
  }  // end of getValues

  template <size_type N>
  inline std::span<real, N> Function::getValues(const size_type o) {
    return std::span<real, N>(this->values.data() + this->getDataOffset(o),
                              this->getNumberOfComponents());
  }  // end of getValues

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_FUNCTION_IXX */
