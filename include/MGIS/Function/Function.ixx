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

  template <size_type data_size, size_type data_offset, size_type data_stride>
  inline size_type
  FunctionDataLayout<data_size, data_offset, data_stride>::getDataOffset(
      const size_type i) const noexcept {
    if constexpr ((data_stride != dynamic_extent) &&
                  (data_offset != dynamic_extent)) {
      if constexpr (data_offset == 0) {
        return i * data_stride;
      } else {
        return i * data_stride + data_offset;
      }
    } else if constexpr (data_offset != dynamic_extent) {
      if constexpr (data_offset == 0) {
        return i * (this->getDataStride());
      } else {
        return i * (this->getDataStride()) + data_offset;
      }
    } else {
      return i * (this->getDataStride()) + this->getDataOffset();
    }
  }  // end of getDataOffset

  template <size_type data_size, size_type data_offset, size_type data_stride>
  constexpr bool has_dynamic_properties(
      const FunctionDataLayout<data_size, data_offset, data_stride>&){
    return (data_size == dynamic_extent) ||    //
           (data_offset == dynamic_extent) ||  //
           (data_stride == dynamic_extent);
  } // end of has_dynamic_properties

  inline const AbstractSpace& ImmutableFunctionView::getSpace() const noexcept {
    return *(this->qspace);
  }  // end of getSpace

  inline std::shared_ptr<const AbstractSpace>
  ImmutableFunctionView::getSpacePointer() const noexcept {
    return this->qspace;
  }  // end of getSpacePointer

  inline const real& ImmutableFunctionView::getIntegrationPointValue(
      const size_type o) const {
    return *(this->immutable_values.data() + this->getDataOffset(o));
  }  // end of getIntegrationPointValues

  inline std::span<const real> ImmutableFunctionView::getIntegrationPointValues(
      const size_type o) const {
    return std::span<const real>(
        this->immutable_values.data() + this->getDataOffset(o),
        this->data_size);
  }  // end of getIntegrationPointValues

  template <size_type N>
  inline std::span<const real, N>
  ImmutableFunctionView::getIntegrationPointValues(const size_type o) const {
    return std::span<const real, N>(
        this->immutable_values.data() + this->getDataOffset(o),
        this->data_size);
  }  // end of getIntegrationPointValues

  inline std::span<const real> ImmutableFunctionView::getValues() const {
    return this->immutable_values;
  }

  inline real& Function::getIntegrationPointValue(const size_type o) {
    return *(this->values.data() + this->getDataOffset(o));
  }  // end of getIntegrationPointValues

  inline std::span<real> Function::getIntegrationPointValues(
      const size_type o) {
    return std::span<real>(this->values.data() + this->getDataOffset(o),
                           this->data_size);
  }  // end of getIntegrationPointValues

  template <size_type N>
  inline std::span<real, N> Function::getIntegrationPointValues(
      const size_type o) {
    return std::span<real, N>(this->values.data() + this->getDataOffset(o),
                              this->data_size);
  }  // end of getIntegrationPointValues

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_FUNCTION_IXX */
