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
    return this->data_stride;
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

  template <std::integral IntegerType>
  void check_positivity(
      IntegerType v,
      const char* const error_message) requires(std::is_signed_v<IntegerType>) {
    raise_if(v < 0, error_message);
  }

  template <std::integral IntegerType>
  void check_positivity(IntegerType, const char* const) noexcept
      requires(std::is_unsigned_v<IntegerType>) {}

  template <FunctionalSpaceConcept Space, DataLayoutDescription layout>
  ImmutableFunctionView<Space, layout>::ImmutableFunctionView() = default;

  template <FunctionalSpaceConcept Space, DataLayoutDescription layout>
  ImmutableFunctionView<Space, layout>::ImmutableFunctionView(
      std::shared_ptr<const Space> s,
      const size_type nv,
      const size_type ds) requires((layout.size != dynamic_extent) &&
                                   (layout.stride != dynamic_extent))
      : space(s) {
    this->data_stride = ds;
    this->data_size = nv;
    raise_if(this->data_size <= 0,
             "ImmutableFunctionView::"
             "ImmutableFunctionView: "
             "invalid data size");
    raise_if(this->data_size > this->data_stride,
             "ImmutableFunctionView::"
             "ImmutableFunctionView: "
             "invalid data range is outside the stride size");
  }  // end of ImmutableFunctionView

  template <FunctionalSpaceConcept Space, DataLayoutDescription layout>
  ImmutableFunctionView<Space, layout>::ImmutableFunctionView(
      std::shared_ptr<const Space> s,
      std::span<const real> v,
      const size_type ds) requires(layout.size == dynamic_extent)
      : space(s) {
    this->data_size = ds;
    raise_if(this->data_size <= 0,
             "Function::Function: invalid "
             "data size");
    if (this->space->size() == 0) {
      // this may happen due to partionning in parallel
      if constexpr (layout.stride != dynamic_extent) {
        this->data_stride = dynamic_extent;
      }
      return;
    }
    const auto quotient =
        static_cast<size_type>(v.size()) / this->space->size();
    const auto remainder =
        static_cast<size_type>(v.size()) % this->space->size();
    raise_if((remainder != 0) || (quotient <= 0),
             "Function::Function: invalid "
             "values size");
    if constexpr (layout.stride != dynamic_extent) {
      this->data_stride = quotient;
      if (this->data_size == dynamic_extent) {
        this->data_size = this->data_stride;
      } else {
        raise_if(this->data_size > this->data_stride,
                 "ImmutableFunctionView::ImmutableFunctionView: invalid "
                 "data range is outside the stride size");
      }
    } else {
    }
    this->immutable_values = v;
  }  // end of ImmutableFunctionView

  template <FunctionalSpaceConcept Space, DataLayoutDescription layout>
  ImmutableFunctionView<Space, layout>::ImmutableFunctionView(
      std::shared_ptr<const Space> s, std::span<const real> v)
      : space(s), immutable_values(v) {
    if (this->space->size() == 0) {
      // this may happen due to partionning in parallel
      if constexpr (layout.stride != dynamic_extent) {
        this->data_stride = dynamic_extent;
      }
      return;
    }
    const auto quotient =
        static_cast<size_type>(v.size()) / this->space->size();
    const auto remainder =
        static_cast<size_type>(v.size()) % this->space->size();
    raise_if((remainder != 0) || (quotient <= 0),
             "ImmutableFunctionView::ImmutableFunctionView: invalid "
             "values size");
    if constexpr (layout.size == dynamic_extent) {
      this->data_size = quotient;
    } else {
      raise_if(layout.size > quotient,
               "ImmutableFunctionView::ImmutableFunctionView: size"
               "does not match");
    }
    if constexpr (layout.stride == dynamic_extent) {
      this->data_stride = quotient;
    } else {
      raise_if(quotient != layout.stride,
               "ImmutableFunctionView::ImmutableFunctionView: stride "
               "does not match");
    }
  }

  template <FunctionalSpaceConcept Space, DataLayoutDescription layout>
  bool ImmutableFunctionView<Space, layout>::checkCompatibility(
      const ImmutableFunctionView& v) const {
    if (this->getSpacePointer() != v.getSpacePointer()) {
      return false;
    }
    return this->data_size == v.getNumberOfComponents();
  }  // end of checkCompatibility

  template <FunctionalSpaceConcept Space, DataLayoutDescription layout>
  const Space& ImmutableFunctionView<Space, layout>::getSpace() const noexcept {
    return *(this->space);
  }  // end of getSpace

  template <FunctionalSpaceConcept Space, DataLayoutDescription layout>
  std::shared_ptr<const Space>
  ImmutableFunctionView<Space, layout>::getSpacePointer() const noexcept {
    return this->space;
  }  // end of getSpacePointer

  template <FunctionalSpaceConcept Space, DataLayoutDescription layout>
  const real& ImmutableFunctionView<Space, layout>::getValue(
      const size_type o) const requires(LinearSpaceConcept<Space>) {
    return *(this->immutable_values.data() + this->getDataOffset(o));
  }  // end of getValues

  template <FunctionalSpaceConcept Space, DataLayoutDescription layout>
  std::span<const real> ImmutableFunctionView<Space, layout>::getValues(
      const size_type o) const requires(LinearSpaceConcept<Space>) {
    return std::span<const real>(
        this->immutable_values.data() + this->getDataOffset(o),
        this->getNumberOfComponents());
  }  // end of getValues

  template <FunctionalSpaceConcept Space, DataLayoutDescription layout>
  template <size_type N>
  std::span<const real, N> ImmutableFunctionView<Space, layout>::getValues(
      const size_type o) const requires(LinearSpaceConcept<Space>) {
    return std::span<const real, N>(
        this->immutable_values.data() + this->getDataOffset(o),
        this->getNumberOfComponents());
  }  // end of getValues

  template <FunctionalSpaceConcept Space, DataLayoutDescription layout>
  std::span<const real> ImmutableFunctionView<Space, layout>::getValues()
      const {
    return this->immutable_values;
  }

  template <FunctionalSpaceConcept Space, DataLayoutDescription layout>
  ImmutableFunctionView<Space, layout>::~ImmutableFunctionView() = default;

  /*!
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
  */

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_FUNCTION_IXX */
