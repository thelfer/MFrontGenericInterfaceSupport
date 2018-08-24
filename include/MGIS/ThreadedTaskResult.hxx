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

#include <exception>
#include <type_traits>
#include "MGIS/Config.hxx"

namespace mgis {

  //! non-template base class of the ThreadedTaskResult class
  struct MGIS_VISIBILITY_EXPORT ThreadedTaskResultBase {
   protected:
    //! thr
    MGIS_NORETURN static void throwBadCastException();
    //! thr
    MGIS_NORETURN static void throwNullException();
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
    MGIS_INLINE ThreadedTaskResult();
    /*!
     * \brief constructor of T
     * \param[in] args: arguments to T constructor
     */
    template <typename... Args>
    MGIS_INLINE ThreadedTaskResult(Args&&...);
    /*!
     * \brief constructor of T
     */
    MGIS_INLINE ThreadedTaskResult(const T&);
    //! \brief move constructor
    MGIS_INLINE ThreadedTaskResult(ThreadedTaskResult&&);
    //! \brief copy constructor
    MGIS_INLINE ThreadedTaskResult(const ThreadedTaskResult&);
    //! \brief assignement
    MGIS_INLINE ThreadedTaskResult& operator=(const ThreadedTaskResult&);
    //! \brief move assignement
    MGIS_INLINE ThreadedTaskResult& operator=(ThreadedTaskResult&&);
    //! \brief assignement
    MGIS_INLINE ThreadedTaskResult& operator=(const T&);
    //! \brief move assignement
    MGIS_INLINE ThreadedTaskResult& operator=(T&&);
    //! \brief set current exception
    MGIS_INLINE void setException(const std::exception_ptr&);
//! \brief throw the catched exception
#ifndef _MSC_VER
    MGIS_NORETURN MGIS_INLINE void rethrow();
#else  /* _MSC_VER */
    MGIS_INLINE void rethrow();
#endif /* _MSC_VER */
       //! \brief conversion to bool
    MGIS_INLINE operator bool() const;
    //! \brief conversion to underlying type
    MGIS_INLINE T& operator*();
    //! \brief conversion to underlying type
    MGIS_INLINE const T& operator*() const;
    //! \brief conversion to underlying type
    MGIS_INLINE T* operator->();
    //! \brief conversion to underlying type
    MGIS_INLINE const T* operator->() const;
    //! \brief destructor
    MGIS_INLINE ~ThreadedTaskResult();

   private:
    /*!
     * build an object
     * \param[in] args: arguments to the constructor
     */
    template <typename... Args>
    MGIS_INLINE void build(Args&&...);
    //! \brief clear the underlying object destructor
    MGIS_INLINE void clear();
    //! \brief pointer the underlying object storage
    MGIS_INLINE T* get_pointer();
    //! \brief pointer the underlying object storage
    MGIS_INLINE const T* get_pointer() const;
    //! \brief if true, the ThreadedTaskResult contains a value
    bool initialized = false;
    //! \brief storage of the underlying value
    typename std::aligned_storage<sizeof(T), alignof(T)>::type storage;
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
    MGIS_NORETURN void rethrow();
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
