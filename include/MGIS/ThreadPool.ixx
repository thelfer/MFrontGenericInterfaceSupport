/*!
 * \file   ThreadPool.ixx
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

#ifndef MGIS_THREAD_POOL_IXX
#define MGIS_THREAD_POOL_IXX

#include <memory>
#include <type_traits>

namespace mgis {

  template <typename F>
  struct ThreadPool::Wrapper {
    Wrapper(F&& f_) : f(f_) {}
    template <typename... Args>
    ThreadedTaskResult<typename std::result_of<F(Args...)>::type> operator()(
        Args&&... args) {
      using result = typename std::result_of<F(Args...)>::type;
      using apply = typename std::conditional<std::is_same<result, void>::value,
                                              GetVoid, Get<result>>::type;
      ThreadedTaskResult<result> r;
      apply::exe(r, f, std::forward<Args>(args)...);
      return r;
    }

   private:
    template <typename T, typename... Args>
    struct Get {
      static void exe(ThreadedTaskResult<T>& r, F& t, Args&&... args) {
        try {
          r = t(args...);
        } catch (...) {
          r.setException(std::current_exception());
        }
      }
    };
    struct GetVoid {
      template <typename... Args>
      static void exe(ThreadedTaskResult<void>& r, F& t, Args&&... args) {
        try {
          t(args...);
        } catch (...) {
          r.setException(std::current_exception());
        }
      }
    };
    // task to be performed
    F f;
  };

  // add new work item to the pool
  template <typename F, typename... Args>
  std::future<ThreadedTaskResult<typename std::result_of<F(Args...)>::type>>
  ThreadPool::addTask(F&& f, Args&&... a) {
    using return_type =
        ThreadedTaskResult<typename std::result_of<F(Args...)>::type>;
    using task = std::packaged_task<return_type()>;
    auto t = std::make_shared<task>(
        std::bind(Wrapper<F>(std::forward<F>(f)), std::forward<Args>(a)...));
    auto res = t->get_future();
    {
      std::unique_lock<std::mutex> lock(this->m);
      // don't allow enqueueing after stopping the pool
      if (this->stop) {
        throw std::runtime_error(
            "ThreadPool::addTask: "
            "enqueue on stopped ThreadPool");
      }
      this->tasks.emplace([t] { (*t)(); });
    }
    this->c.notify_one();
    return res;
  }

}  // end of namespace mgis

#endif /* MGIS_THREAD_POOL_IXX */
