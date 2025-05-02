/*!
 * \file   MGIS/QuadratureFunction/QuadratureFunction.ixx
 * \brief
 * \author Thomas Helfer
 * \date   10/03/2021
 */

#ifndef LIB_MGIS_QUADRATUREFUNCTION_QUADRATUREFUNCTION_IXX
#define LIB_MGIS_QUADRATUREFUNCTION_QUADRATUREFUNCTION_IXX

namespace mgis::quadrature_function {

  template <size_type offset>
  constexpr size_type QuadratureFunctionDataOffset<offset>::getDataOffset()
      const noexcept {
    return offset;
  }  // end of getDataOffset

  inline size_type QuadratureFunctionDataOffset<dynamic_extent>::getDataOffset()
      const noexcept {
    return this->data_begin;
  }  // end of getDataOffset

  template <size_type data_size>
  constexpr bool QuadratureFunctionDataSize<data_size>::isScalar() const noexcept {
    return data_size == 1;
  }  // end of isScalar

  template <size_type data_size>
  constexpr size_type QuadratureFunctionDataSize<data_size>::getNumberOfComponents()
      const noexcept {
    return data_size;
  }  // end of QuadratureFunctionDataSize<dynamic_extent>::getNumberOfComponents

  inline bool QuadratureFunctionDataSize<dynamic_extent>::isScalar() const noexcept {
    return this->getNumberOfComponents() == 1;
  }  // end of isScalar

  inline size_type QuadratureFunctionDataSize<dynamic_extent>::getNumberOfComponents()
      const noexcept {
    return this->data_size;
  }  // end of QuadratureFunctionDataSize<dynamic_extent>::getNumberOfComponents

  template <size_type data_stride>
  constexpr size_type QuadratureFunctionDataStride<data_stride>::getDataStride()
      const noexcept {
    return this->data_stride;
  }  // end of getDataStride

  inline size_type QuadratureFunctionDataStride<dynamic_extent>::getDataStride()
      const noexcept {
    return this->data_stride;
  }  // end of getDataStride

  inline size_type QuadratureFunctionDataLayout::getDataOffset(
      const size_type i) const noexcept {
    return i * (this->data_stride) + this->getDataOffset();
  }  // end of getDataOffset

  inline const AbstractQuadratureSpace&
  ImmutableQuadratureFunctionView::getQuadratureSpace() const noexcept {
    return *(this->qspace);
  }  // end of getQuadratureSpace

  inline std::shared_ptr<const AbstractQuadratureSpace>
  ImmutableQuadratureFunctionView::getQuadratureSpacePointer() const noexcept {
    return this->qspace;
  }  // end of getQuadratureSpacePointer

  inline const real&
  ImmutableQuadratureFunctionView::getIntegrationPointValue(
      const size_type o) const {
    return *(this->immutable_values.data() + this->getDataOffset(o));
  }  // end of getIntegrationPointValues

  inline std::span<const real>
  ImmutableQuadratureFunctionView::getIntegrationPointValues(
      const size_type o) const {
    return std::span<const real>(
        this->immutable_values.data() + this->getDataOffset(o),
        this->data_size);
  }  // end of getIntegrationPointValues

  template <size_type N>
  inline std::span<const real, N>
  ImmutableQuadratureFunctionView::getIntegrationPointValues(
      const size_type o) const {
    return std::span<const real, N>(
        this->immutable_values.data() + this->getDataOffset(o),
        this->data_size);
  }  // end of getIntegrationPointValues

  inline std::span<const real>
  ImmutableQuadratureFunctionView::getValues() const {
    return this->immutable_values;
  }
  
  inline real& QuadratureFunction::getIntegrationPointValue(
      const size_type o) {
    return *(this->values.data() + this->getDataOffset(o));
  }  // end of getIntegrationPointValues

  inline std::span<real> QuadratureFunction::getIntegrationPointValues(
      const size_type o) {
    return std::span<real>(this->values.data() + this->getDataOffset(o),
                            this->data_size);
  }  // end of getIntegrationPointValues

  template <size_type N>
  inline std::span<real, N>
  QuadratureFunction::getIntegrationPointValues(const size_type o) {
    return std::span<real, N>(this->values.data() + this->getDataOffset(o),
                               this->data_size);
  }  // end of getIntegrationPointValues

}  // end of namespace mgis::quadrature_function

#endif /* LIB_MGIS_QUADRATUREFUNCTION_QUADRATUREFUNCTION_IXX */
