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
  ThreadedTaskResult<T>::ThreadedTaskResult(Args&&... args) {
    this->build(args...);
  }  // end of ThreadedTaskResult<T>::ThreadedTaskResult

  template <typename T>
  ThreadedTaskResult<T>::ThreadedTaskResult(ThreadedTaskResult&& src) {
    this->eptr = std::move(src.eptr);
    if (src.initialized) {
      this->initialized = true;
      this->build(std::move(*(src.get_pointer())));
    }
  }

  template <typename T>
  ThreadedTaskResult<T>::ThreadedTaskResult(const ThreadedTaskResult& src) {
    this->eptr = src.eptr;
    if (src.initialized) {
      this->initialized = true;
      this->build(src.get_pointer());
    }
  }

  template <typename T>
  ThreadedTaskResult<T>& ThreadedTaskResult<T>::operator=(
      ThreadedTaskResult&& src) {
    if (this != &src) {
      this->eptr = std::move(src.eptr);
      if (src.initialized) {
        if (this->initialized) {
          *(this->get_pointer()) = std::move(*(src.get_pointer()));
        } else {
          this->build(std::move(*(src.get_pointer())));
        }
      } else {
        this->clear();
      }
    }
    return *this;
  }

  template <typename T>
  ThreadedTaskResult<T>& ThreadedTaskResult<T>::operator=(
      const ThreadedTaskResult& src) {
    if (this != &src) {
      this->eptr = src.eptr;
      if (src.initialized) {
        if (this->initialized) {
          *(this->get_pointer()) = *(src.get_pointer());
        } else {
          this->build(*(src.get_pointer()));
        }
      } else {
        this->clear();
      }
    }
    return *this;
  }

  template <typename T>
  ThreadedTaskResult<T>& ThreadedTaskResult<T>::operator=(const T& src) {
    if (this->initialized) {
      *(this->get_pointer()) = src;
    } else {
      this->build(src);
    }
    return *this;
  }

  template <typename T>
  ThreadedTaskResult<T>& ThreadedTaskResult<T>::operator=(T&& src) {
    if (this->initialized) {
      *(this->get_pointer()) = std::forward<T>(src);
    } else {
      this->build(std::forward<T>(src));
    }
    return *this;
  }

  template <typename T>
  ThreadedTaskResult<T>::operator bool() const {
    return this->initialized && (this->eptr == nullptr);
  }

  template <typename T>
  T& ThreadedTaskResult<T>::operator*() {
    if (this->eptr != nullptr) {
      this->rethrow();
    }
    if (!this->initialized) {
      ThreadedTaskResultBase::throwBadCastException();
    }
    return *(this->get_pointer());
  }

  template <typename T>
  const T& ThreadedTaskResult<T>::operator*() const {
    if (this->eptr != nullptr) {
      this->rethrow();
    }
    if (!this->initialized) {
      ThreadedTaskResultBase::throwBadCastException();
    }
    return *(this->get_pointer());
  }

  template <typename T>
  T* ThreadedTaskResult<T>::operator->() {
    if (this->eptr != nullptr) {
      this->rethrow();
    }
    if (!this->initialized) {
      ThreadedTaskResultBase::throwBadCastException();
    }
    return this->get_pointer();
  }

  template <typename T>
  const T* ThreadedTaskResult<T>::operator->() const {
    if (this->eptr != nullptr) {
      this->rethrow();
    }
    if (!this->initialized) {
      ThreadedTaskResultBase::throwBadCastException();
    }
    return this->get_pointer();
  }

  template <typename T>
  void ThreadedTaskResult<T>::setException(const std::exception_ptr& e) {
    this->clear();
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
  template <typename... Args>
  void ThreadedTaskResult<T>::build(Args&&... args) {
    this->clear();
    new (this->get_pointer()) T(args...);
    this->initialized = true;
  }

  template <typename T>
  void ThreadedTaskResult<T>::clear() {
    this->eptr = nullptr;
    if (this->initialized) {
      auto ptr = this->get_pointer();
      ptr->~T();
      this->initialized = false;
    }
  }  // end of clear

  template <typename T>
  T* ThreadedTaskResult<T>::get_pointer() {
    return reinterpret_cast<T*>(::std::addressof(this->storage));
  }  // end of get_pointer

  template <typename T>
  const T* ThreadedTaskResult<T>::get_pointer() const {
    return reinterpret_cast<const T*>(::std::addressof(this->storage));
  }  // end of get_pointer

  template <typename T>
  ThreadedTaskResult<T>::~ThreadedTaskResult() {
    this->clear();
  }

}  // end of namespace mgis

#endif /* LIB_MGIS_THREADEDTASKRESULT_IXX */
