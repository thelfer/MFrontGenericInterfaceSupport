/*!
 * \file   src/QuadratureFunction.cxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <type_traits>
#include "MGIS/Raise.hxx"
#include "MGIS/QuadratureFunction/AbstractQuadratureSpace.hxx"
#include "MGIS/QuadratureFunction/QuadratureFunction.hxx"

namespace mgis::quadrature_function {

  template <std::integral IntegerType>
  void check_positivity(
      IntegerType v,
      const char* const error_message) requires(std::is_signed_v<IntegerType>) {
    raise_if(v < 0, error_message);
  }

  template <std::integral IntegerType>
  void check_positivity(IntegerType, const char* const) noexcept
      requires(std::is_unsigned_v<IntegerType>) {}

  ImmutableQuadratureFunctionView::ImmutableQuadratureFunctionView() = default;

  ImmutableQuadratureFunctionView::ImmutableQuadratureFunctionView(
      std::shared_ptr<const AbstractQuadratureSpace> s,
      const size_type nv,
      const size_type db,
      const size_type ds)
      : qspace(s) {
    if (this->qspace->getNumberOfIntegrationPoints() == 0) {
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
                     "ImmutableQuadratureFunctionView::"
                     "ImmutableQuadratureFunctionView: "
                     "invalid start of the data");
    raise_if(this->data_size <= 0,
             "ImmutableQuadratureFunctionView::"
             "ImmutableQuadratureFunctionView: "
             "invalid data size");
    raise_if(this->data_begin + this->data_size > this->data_stride,
             "ImmutableQuadratureFunctionView::"
             "ImmutableQuadratureFunctionView: "
             "invalid data range is outside the stride size");
  }  // end of ImmutableQuadratureFunctionView

  ImmutableQuadratureFunctionView::ImmutableQuadratureFunctionView(
      std::shared_ptr<const AbstractQuadratureSpace> s,
      std::span<const real> v,
      const size_type db,
      const size_type ds)
      : qspace(s) {
    if (this->qspace->getNumberOfIntegrationPoints() == 0) {
      // this may happen due to partionning in parallel
      this->data_begin = 0;
      this->data_stride = 0;
      this->data_size = 0;
      return;
    }
    this->data_begin = db;
    this->data_size = ds;
    check_positivity(this->data_begin,
                     "ImmutableQuadratureFunctionView::"
                     "ImmutableQuadratureFunctionView: "
                     "invalid start of the data");
    raise_if(this->data_size <= 0,
             "QuadratureFunction::QuadratureFunction: invalid "
             "data size");
    const auto quotient = static_cast<size_type>(v.size()) /
                          this->qspace->getNumberOfIntegrationPoints();
    const auto remainder = static_cast<size_type>(v.size()) %
                           this->qspace->getNumberOfIntegrationPoints();
    raise_if((remainder != 0) || (quotient <= 0),
             "QuadratureFunction::QuadratureFunction: invalid "
             "values size");
    this->data_stride = quotient;
    raise_if(this->data_begin >= this->data_stride,
             "QuadratureFunction::QuadratureFunction: invalid "
             "start of the data");
    if (ds == std::numeric_limits<size_type>::max()) {
      this->data_size = this->data_stride - this->data_begin;
    } else {
      raise_if(this->data_begin + ds > this->data_stride,
               "QuadratureFunction::QuadratureFunction: invalid "
               "data range is outside the stride size");
      this->data_size = ds;
    }
    this->immutable_values = v;
  }  // end of ImmutableQuadratureFunctionView

  bool ImmutableQuadratureFunctionView::checkCompatibility(
      const ImmutableQuadratureFunctionView& v) const {
    if (this->getQuadratureSpacePointer() != v.getQuadratureSpacePointer()) {
      return false;
    }
    return this->data_size == v.getNumberOfComponents();
  }  // end of checkCompatibility

  ImmutableQuadratureFunctionView::~ImmutableQuadratureFunctionView() = default;

  QuadratureFunction::QuadratureFunction(QuadratureFunction&& f,
                                         const bool local_copy) {
    if (!f.local_values_storage.empty()) {
      // the function holds the memory, just take it from him
      static_cast<QuadratureFunctionDataLayout&>(*this).operator=(f);
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
  }  // end of QuadratureFunction

  QuadratureFunction::QuadratureFunction(const QuadratureFunction& v)
      : ImmutableQuadratureFunctionView() {
    this->copy(v);
  }  // end of QuadratureFunction

  QuadratureFunction::QuadratureFunction(
      const ImmutableQuadratureFunctionView& v) {
    this->copy(v);
  }  // end of ImmutableQuadratureFunctionView

  QuadratureFunction::QuadratureFunction(
      std::shared_ptr<const AbstractQuadratureSpace> s, const size_type nv)
      : ImmutableQuadratureFunctionView(s, nv, 0, nv) {
    this->local_values_storage.resize(
        this->qspace->getNumberOfIntegrationPoints() * this->data_size);
    this->values = std::span<real>(this->local_values_storage);
    this->immutable_values = std::span<const real>(this->local_values_storage);
  }  // end of QuadratureFunction::QuadratureFunction

  QuadratureFunction::QuadratureFunction(
      std::shared_ptr<const AbstractQuadratureSpace> s,
      std::span<real> v,
      const size_type db,
      const size_type ds)
      : ImmutableQuadratureFunctionView(s, v, db, ds),
        values(v) {}  // end of QuadratureFunction::QuadratureFunction

  void QuadratureFunction::makeView(QuadratureFunction& f) {
    static_cast<QuadratureFunctionDataLayout&>(*this).operator=(f);
    this->qspace = f.qspace;
    this->values = f.values;
    this->immutable_values = f.immutable_values;
  }

  QuadratureFunction& QuadratureFunction::operator=(
      const QuadratureFunction& src) {
    if (&src != this) {
      this->copy(src);
    }
    return *this;
  }

  QuadratureFunction& QuadratureFunction::operator=(
      const ImmutableQuadratureFunctionView& src) {
    if (&src != this) {
      this->copy(src);
    }
    return *this;
  }

  void QuadratureFunction::copy(const ImmutableQuadratureFunctionView& v) {
    this->qspace = v.getQuadratureSpacePointer();
    const auto n = this->qspace->getNumberOfIntegrationPoints();
    this->data_begin = size_type{};
    this->data_size = v.getNumberOfComponents();
    this->data_stride = v.getNumberOfComponents();
    this->local_values_storage.resize(this->data_size * n);
    this->values = local_values_storage;
    this->immutable_values = local_values_storage;
    this->copyValues(v);
  }  // end of copy

  void QuadratureFunction::copyValues(
      const ImmutableQuadratureFunctionView& v) {
    const auto n =
      this->getQuadratureSpace().getNumberOfIntegrationPoints();
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

  QuadratureFunction::~QuadratureFunction() = default;

}  // end of namespace mgis::quadrature_function
