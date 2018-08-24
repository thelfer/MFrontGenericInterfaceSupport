/*!
 * \file   ThreadedTaskResult.cxx
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

#include <typeinfo>
#include <stdexcept>
#include "MGIS/Raise.hxx"
#include "MGIS/ThreadedTaskResult.hxx"

namespace mgis {

  void ThreadedTaskResultBase::throwBadCastException() {
    throw(std::bad_cast());
  }  // end of ThreadPool::ResultBase::throwBadCastException

  void ThreadedTaskResultBase::throwNullException() {
    raise<std::runtime_error>(
        "ThreadPool::Result::rethrow: "
        "no exception defined");
  }  // end of ThreadPool::ResultBase::throwNullException

  ThreadedTaskResult<void>::ThreadedTaskResult() = default;
  ThreadedTaskResult<void>::ThreadedTaskResult(ThreadedTaskResult&&) = default;
  ThreadedTaskResult<void>::ThreadedTaskResult(const ThreadedTaskResult&) =
      default;
  ThreadedTaskResult<void>& ThreadedTaskResult<void>::operator=(
      ThreadedTaskResult&&) = default;
  ThreadedTaskResult<void>& ThreadedTaskResult<void>::operator=(
      const ThreadedTaskResult&) = default;
  ThreadedTaskResult<void>::~ThreadedTaskResult() = default;

  ThreadedTaskResult<void>::operator bool() const {
    return this->eptr == nullptr;
  }

  void ThreadedTaskResult<void>::setException(const std::exception_ptr& e) {
    this->eptr = e;
  }  // end of ThreadedTaskResult<void>::setException

  void ThreadedTaskResult<void>::rethrow() {
    if (this->eptr == nullptr) {
      ThreadedTaskResultBase::throwNullException();
    }
    std::rethrow_exception(eptr);
  }  // end of ThreadedTaskResult<void>::setException

}  // end of namespace mgis
