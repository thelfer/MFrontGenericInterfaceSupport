/*!
 * \file   ThreadedTaskResult.ixx
 * \brief
 * \author Thomas Helfer
 * \date   24/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_THREADEDTASKRESULT_IXX
#define LIB_MGIS_THREADEDTASKRESULT_IXX

#include <memory>
#include <utility>

namespace mgis {

  template <typename T>
  ThreadedTaskResult<T>::ThreadedTaskResult() = default;

  template <typename T>
  template <typename... Args>
  ThreadedTaskResult<T>::ThreadedTaskResult(Args&&... args)
      : result(args...) {}  // end of ThreadedTaskResult

  template <typename T>
  ThreadedTaskResult<T>::ThreadedTaskResult(ThreadedTaskResult&&) = default;

  template <typename T>
  ThreadedTaskResult<T>::ThreadedTaskResult(const ThreadedTaskResult&) =
      default;

  template <typename T>
  ThreadedTaskResult<T>& ThreadedTaskResult<T>::operator=(
      ThreadedTaskResult&&) = default;

  template <typename T>
  ThreadedTaskResult<T>& ThreadedTaskResult<T>::operator=(
      const ThreadedTaskResult& src) = default;

  template <typename T>
  ThreadedTaskResult<T>& ThreadedTaskResult<T>::operator=(const T& src) {
    this->result = src;
    return *this;
  }

  template <typename T>
  ThreadedTaskResult<T>& ThreadedTaskResult<T>::operator=(T&& src) {
    this->result = std::forward<T>(src);
    return *this;
  }

  template <typename T>
  ThreadedTaskResult<T>::operator bool() const {
    return this->result.has_value() && (this->eptr == nullptr);
  }

  template <typename T>
  T& ThreadedTaskResult<T>::operator*() {
    if (this->eptr != nullptr) {
      this->rethrow();
    }
    if (!this->result.has_value()) {
      ThreadedTaskResultBase::throwBadCastException();
    }
    return *(this->result);
  }  // end of operator *

  template <typename T>
  const T& ThreadedTaskResult<T>::operator*() const {
    if (this->eptr != nullptr) {
      this->rethrow();
    }
    if (!this->result_has_value()) {
      ThreadedTaskResultBase::throwBadCastException();
    }
    return *(this->result);
  }

  template <typename T>
  T* ThreadedTaskResult<T>::operator->() {
    if (this->eptr != nullptr) {
      this->rethrow();
    }
    if (!this->result.has_value()) {
      ThreadedTaskResultBase::throwBadCastException();
    }
    return this->result.operator->();
  }

  template <typename T>
  const T* ThreadedTaskResult<T>::operator->() const {
    if (this->eptr != nullptr) {
      this->rethrow();
    }
    if (!this->result.has_value()) {
      ThreadedTaskResultBase::throwBadCastException();
    }
    return this->result.operator->();
  }

  template <typename T>
  void ThreadedTaskResult<T>::setException(const std::exception_ptr& e) {
    this->result.reset();
    this->eptr = e;
  }  // end of ThreadedTaskResult<T>::setException

  template <typename T>
  void ThreadedTaskResult<T>::rethrow() {
    if (this->eptr == nullptr) {
      ThreadedTaskResultBase::throwNullException();
    }
    std::rethrow_exception(eptr);
  }  // end of ThreadedTaskResult<T>::setException

  template <typename T>
  ThreadedTaskResult<T>::~ThreadedTaskResult() = default;

}  // end of namespace mgis

#endif /* LIB_MGIS_THREADEDTASKRESULT_IXX */
