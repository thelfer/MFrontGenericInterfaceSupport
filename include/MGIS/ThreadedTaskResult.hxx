/*!
 * \file   ThreadedTaskResult.hxx
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

#ifndef LIB_MGIS_THREADEDTASKRESULT_HXX
#define LIB_MGIS_THREADEDTASKRESULT_HXX

#include <optional>
#include <exception>
#include <type_traits>
#include "MGIS/Config.hxx"

namespace mgis {

  //! non-template base class of the ThreadedTaskResult class
  struct MGIS_VISIBILITY_EXPORT ThreadedTaskResultBase {
   protected:
    //! thr
    [[noreturn]] static void throwBadCastException();
    //! thr
    [[noreturn]] static void throwNullException();
  };  // end of struct ThreadedTaskResultBase

  /*!
   * \brief a class standing for the result of a taks
   *
   * Its interface is loosely base on the interface of
   * std::optional (not available in C++11), but it also stores an
   * std::exception_ptr if .
   */
  template <typename T>
  struct ThreadedTaskResult : public ThreadedTaskResultBase {
    //! \brief default constructor
    inline ThreadedTaskResult();
    /*!
     * \brief constructor of T
     * \param[in] args: arguments to T constructor
     */
    template <typename... Args>
    inline ThreadedTaskResult(Args&&...);
    /*!
     * \brief constructor of T
     */
    inline ThreadedTaskResult(const T&);
    //! \brief move constructor
    inline ThreadedTaskResult(ThreadedTaskResult&&);
    //! \brief copy constructor
    inline ThreadedTaskResult(const ThreadedTaskResult&);
    //! \brief assignement
    inline ThreadedTaskResult& operator=(const ThreadedTaskResult&);
    //! \brief move assignement
    inline ThreadedTaskResult& operator=(ThreadedTaskResult&&);
    //! \brief assignement
    inline ThreadedTaskResult& operator=(const T&);
    //! \brief move assignement
    inline ThreadedTaskResult& operator=(T&&);
    //! \brief set current exception
    inline void setException(const std::exception_ptr&);
//! \brief throw the catched exception
#ifndef _MSC_VER
    [[noreturn]] inline void rethrow();
#else  /* _MSC_VER */
    inline void rethrow();
#endif /* _MSC_VER */
       //! \brief conversion to bool
    inline operator bool() const;
    //! \brief conversion to underlying type
    inline T& operator*();
    //! \brief conversion to underlying type
    inline const T& operator*() const;
    //! \brief conversion to underlying type
    inline T* operator->();
    //! \brief conversion to underlying type
    inline const T* operator->() const;
    //! \brief destructor
    inline ~ThreadedTaskResult();

   private:
    //! \brief result of the ThreadedTaskResult
    std::optional<T> result;
    //! exception ptr thrown during the task
    std::exception_ptr eptr;
  };

  //! Partial specialisation for non-returning tasks
  template <>
  struct MGIS_VISIBILITY_EXPORT ThreadedTaskResult<void>
      : public ThreadedTaskResultBase {
    //! \brief default constructor
    ThreadedTaskResult();
    //! \brief move constructor
    ThreadedTaskResult(ThreadedTaskResult&&);
    //! \brief copy constructor
    ThreadedTaskResult(const ThreadedTaskResult&);
    //! \brief assignement
    ThreadedTaskResult& operator=(const ThreadedTaskResult&);
    //! \brief move assignement
    ThreadedTaskResult& operator=(ThreadedTaskResult&&);
    //! \brief set current exception
    void setException(const std::exception_ptr&);
    //! \brief throw the catched exception
    [[noreturn]] void rethrow();
    //! \brief conversion to bool
    operator bool() const;
    //! \brief destructor
    ~ThreadedTaskResult();

   private:
    //! exception ptr thrown during the task
    std::exception_ptr eptr;
  };

}  // end of namespace mgis

#include "MGIS/ThreadedTaskResult.ixx"

#endif /* LIB_MGIS_THREADEDTASKRESULT_HXX */
