/*!
 * \file   ThreadPool.cxx
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

#include <memory>
#include <stdexcept>
#include "MGIS/ThreadPool.hxx"

namespace mgis {

  ThreadPool::ThreadPool(const size_t n) {
    this->statuses.resize(n, ThreadPool::Status::IDLE);
    for (size_t i = 0; i < n; ++i) {
      auto f = [this, i] {
        for (;;) {
          std::function<void()> task;
          {
            std::unique_lock<std::mutex> lock(this->m);
            while (!(this->stop || !this->tasks.empty())) {
              this->c.wait(
                  lock, [this] { return this->stop || !this->tasks.empty(); });
            }
            if (this->stop && this->tasks.empty()) {
              return;
            }
            task = std::move(this->tasks.front());
            this->tasks.pop();
            this->statuses[i] = ThreadPool::Status::WORKING;
            this->c.notify_all();
          }
          task();
          {
            std::unique_lock<std::mutex> lock(this->m);
            this->statuses[i] = ThreadPool::Status::IDLE;
            this->c.notify_all();
          }
        }
      };
      this->workers.emplace_back(f);
    }
  }

  size_type ThreadPool::getNumberOfThreads() const {
    return this->workers.size();
  }  // end of ThreadPool::getNumberOfThreads

  void ThreadPool::wait() {
    std::unique_lock<std::mutex> lock(this->m);
    while (!this->tasks.empty()) {
      this->c.wait(lock, [this] { return this->tasks.empty(); });
    }
    for (decltype(this->statuses.size()) i = 0; i != this->statuses.size();
         ++i) {
      while (this->statuses[i] != Status::IDLE) {
        this->c.wait(lock,
                     [this, i] { return this->statuses[i] == Status::IDLE; });
      }
    }
  }  // end of ThreadPool::wait()

  ThreadPool::~ThreadPool() {
    {
      std::unique_lock<std::mutex> lock(m);
      this->stop = true;
    }
    this->c.notify_all();
    // the destructor joins all threads
    for (auto &w : this->workers) {
      w.join();
    }
  }

}  // end of namespace mgis
