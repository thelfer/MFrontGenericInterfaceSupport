/*!
 * \file   ThreadPool.hxx
 * \brief A ThreadPool implementation based on the initial
 * implementation of Jakob Progsch, Václav Zeman:
 * <https://github.com/progschj/ThreadPool>
 *
 * We added the possibility to handle exceptions through the
 * ThreadedTaskResult class.
 *
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

#ifndef MGIS_THREAD_POOL_HXX
#define MGIS_THREAD_POOL_HXX

#include <queue>
#include <mutex>
#include <vector>
#include <thread>
#include <future>
#include <functional>
#include <condition_variable>
#include "MGIS/Config.hxx"
#include "MGIS/ThreadedTaskResult.hxx"

namespace mgis {

  /*!
   * \brief structure handling a fixed-size pool of threads
   */
  struct MGIS_VISIBILITY_EXPORT ThreadPool {
    /*!
     * \brief constructor
     * \param[in] n: number of thread to be created
     */
    ThreadPool(const size_type);
    /*!
     * \brief add a new task
     * \param[in] f: task
     * \param[in] a: arguments passed to the the task
     */
    template <typename F, typename... Args>
    std::future<ThreadedTaskResult<typename std::result_of<F(Args...)>::type>>
    addTask(F&&, Args&&...);
    //! \return the number of threads managed by the ppol
    size_type getNumberOfThreads() const;
    //! \brief wait for all tasks to be finished
    void wait();
    //! destructor
    ~ThreadPool();

   private:
    //! wrapper around the given task
    template <typename F>
    struct Wrapper;
    enum Status { WORKING, IDLE };  // end of enum Status
    std::vector<Status> statuses;
    //! list of threads
    std::vector<std::thread> workers;
    // the task queue
    std::queue<std::function<void()>> tasks;
    // synchronization
    std::mutex m;
    std::condition_variable c;
    bool stop = false;
  };

}  // end of namespace mgis

#include "MGIS/ThreadPool.ixx"

#endif /* MGIS_THREAD_POOL_HXX */
