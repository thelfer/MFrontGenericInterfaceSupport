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
#include "MGIS/Function/Space.hxx"
#include "MGIS/Function/Function.hxx"

namespace mgis::function {

  /*
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

    Function::Function(std::shared_ptr<const AbstractSpace> s, const size_type
    nv) : ImmutableFunctionView(s, nv, nv) {
      this->local_values_storage.resize(this->qspace->getSpaceSize() *
                                        this->data_size);
      this->values = std::span<real>(this->local_values_storage);
      this->immutable_values = std::span<const
    real>(this->local_values_storage); }  // end of Function::Function

    Function::Function(std::shared_ptr<const AbstractSpace> s,
                       std::span<real> v,
                       const size_type ds)
        : ImmutableFunctionView(s, v, ds),
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
      const auto* const v_values = v.getValues().data();
      const auto vs = v.getDataStride();
      if (vs == v.getNumberOfComponents()) {
        // data are continous in v
        std::copy(v_values, v_values + this->values.size(),
    this->values.begin()); } else { if (this->data_size == 1) {
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
  */

}  // end of namespace mgis::function
