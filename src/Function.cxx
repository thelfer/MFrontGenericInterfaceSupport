/*!
 * \file   src/Function.cxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <type_traits>
#include "MGIS/Raise.hxx"
#include "MGIS/Function/AbstractSpace.hxx"
#include "MGIS/Function/Function.hxx"

namespace mgis::function {

  template <std::integral IntegerType>
  void check_positivity(
      IntegerType v,
      const char* const error_message) requires(std::is_signed_v<IntegerType>) {
    raise_if(v < 0, error_message);
  }

  template <std::integral IntegerType>
  void check_positivity(IntegerType, const char* const) noexcept
      requires(std::is_unsigned_v<IntegerType>) {}

  ImmutableFunctionView::ImmutableFunctionView() = default;

  ImmutableFunctionView::ImmutableFunctionView(
      std::shared_ptr<const AbstractSpace> s,
      const size_type nv,
      const size_type db,
      const size_type ds)
      : qspace(s) {
    if (this->qspace->getSpaceSize() == 0) {
      // this may happen due to partionning in parallel
      this->data_begin = 0;
      this->data_stride = 0;
      this->data_size = 0;
      return;
    }
    this->data_stride = ds;
    this->data_begin = db;
    this->data_size = nv;
    check_positivity(this->data_begin,
                     "ImmutableFunctionView::"
                     "ImmutableFunctionView: "
                     "invalid start of the data");
    raise_if(this->data_size <= 0,
             "ImmutableFunctionView::"
             "ImmutableFunctionView: "
             "invalid data size");
    raise_if(this->data_begin + this->data_size > this->data_stride,
             "ImmutableFunctionView::"
             "ImmutableFunctionView: "
             "invalid data range is outside the stride size");
  }  // end of ImmutableFunctionView

  ImmutableFunctionView::ImmutableFunctionView(
      std::shared_ptr<const AbstractSpace> s,
      std::span<const real> v,
      const size_type db,
      const size_type ds)
      : qspace(s) {
    if (this->qspace->getSpaceSize() == 0) {
      // this may happen due to partionning in parallel
      this->data_begin = 0;
      this->data_stride = 0;
      this->data_size = 0;
      return;
    }
    this->data_begin = db;
    this->data_size = ds;
    check_positivity(this->data_begin,
                     "ImmutableFunctionView::"
                     "ImmutableFunctionView: "
                     "invalid start of the data");
    raise_if(this->data_size <= 0,
             "Function::Function: invalid "
             "data size");
    const auto quotient =
        static_cast<size_type>(v.size()) / this->qspace->getSpaceSize();
    const auto remainder =
        static_cast<size_type>(v.size()) % this->qspace->getSpaceSize();
    raise_if((remainder != 0) || (quotient <= 0),
             "Function::Function: invalid "
             "values size");
    this->data_stride = quotient;
    raise_if(this->data_begin >= this->data_stride,
             "Function::Function: invalid "
             "start of the data");
    if (ds == std::numeric_limits<size_type>::max()) {
      this->data_size = this->data_stride - this->data_begin;
    } else {
      raise_if(this->data_begin + ds > this->data_stride,
               "Function::Function: invalid "
               "data range is outside the stride size");
      this->data_size = ds;
    }
    this->immutable_values = v;
  }  // end of ImmutableFunctionView

  bool ImmutableFunctionView::checkCompatibility(
      const ImmutableFunctionView& v) const {
    if (this->getSpacePointer() != v.getSpacePointer()) {
      return false;
    }
    return this->data_size == v.getNumberOfComponents();
  }  // end of checkCompatibility

  ImmutableFunctionView::~ImmutableFunctionView() = default;

  Function::Function(Function&& f, const bool local_copy) {
    if (!f.local_values_storage.empty()) {
      // the function holds the memory, just take it from him
      static_cast<DataLayout&>(*this).operator=(f);
      this->qspace = f.qspace;
      this->local_values_storage = std::move(f.local_values_storage);
      this->values = local_values_storage;
      this->immutable_values = local_values_storage;
    } else {
      // the function does not hold the memory
      if (local_copy) {
        this->copy(f);
      } else {
        this->makeView(f);
      }
    }
  }  // end of Function

  Function::Function(const Function& v) : ImmutableFunctionView() {
    this->copy(v);
  }  // end of Function

  Function::Function(const ImmutableFunctionView& v) {
    this->copy(v);
  }  // end of ImmutableFunctionView

  Function::Function(std::shared_ptr<const AbstractSpace> s, const size_type nv)
      : ImmutableFunctionView(s, nv, 0, nv) {
    this->local_values_storage.resize(this->qspace->getSpaceSize() *
                                      this->data_size);
    this->values = std::span<real>(this->local_values_storage);
    this->immutable_values = std::span<const real>(this->local_values_storage);
  }  // end of Function::Function

  Function::Function(std::shared_ptr<const AbstractSpace> s,
                     std::span<real> v,
                     const size_type db,
                     const size_type ds)
      : ImmutableFunctionView(s, v, db, ds),
        values(v) {}  // end of Function::Function

  void Function::makeView(Function& f) {
    static_cast<DataLayout&>(*this).operator=(f);
    this->qspace = f.qspace;
    this->values = f.values;
    this->immutable_values = f.immutable_values;
  }

  Function& Function::operator=(const Function& src) {
    if (&src != this) {
      this->copy(src);
    }
    return *this;
  }

  Function& Function::operator=(const ImmutableFunctionView& src) {
    if (&src != this) {
      this->copy(src);
    }
    return *this;
  }

  void Function::copy(const ImmutableFunctionView& v) {
    this->qspace = v.getSpacePointer();
    const auto n = this->qspace->getSpaceSize();
    this->data_begin = size_type{};
    this->data_size = v.getNumberOfComponents();
    this->data_stride = v.getNumberOfComponents();
    this->local_values_storage.resize(this->data_size * n);
    this->values = local_values_storage;
    this->immutable_values = local_values_storage;
    this->copyValues(v);
  }  // end of copy

  void Function::copyValues(const ImmutableFunctionView& v) {
    const auto n = this->getSpace().getSpaceSize();
    if (n == 0) {
      return;
    }
    const auto* const v_values = v.getValues().data() + v.getDataOffset();
    const auto vs = v.getDataStride();
    if (vs == v.getNumberOfComponents()) {
      // data are continous in v
      std::copy(v_values, v_values + this->values.size(), this->values.begin());
    } else {
      if (this->data_size == 1) {
        // special case for scalars
        for (size_type i = 0; i != this->values.size(); ++i) {
          this->values[i] = v_values[i * vs];
        }
      } else {
        auto pv = this->values.begin();
        for (size_type i = 0; i != n; ++i) {
          const auto b = v_values + i * vs;
          const auto e = b + this->data_size;
          std::copy(b, e, pv);
          pv += this->data_size;
        }
      }
    }
  }  // end of copy

  Function::~Function() = default;

}  // end of namespace mgis::function
