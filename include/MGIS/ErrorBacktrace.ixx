/*!
 * \file   MGIS/ErrorBacktrace.ixx
 * \brief  This file implements the inline methods of the `ErrorBacktrace`
 * class. \date   01/09/2026 \copyright (C) Copyright Thomas Helfer 2018. Use,
 * modification and distribution are subject to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_ERRORBACKTRACE_IXX
#define LIB_MGIS_ERRORBACKTRACE_IXX 1

namespace mgis {

#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION

  inline void ErrorBacktrace::assertOrTerminate(const bool b,
                                                const char *msg,
                                                const std::source_location &l) {
    if (!b) {
      this->terminate(msg, l);
    }
  }  // end of assertOrTerminate

  inline void ErrorBacktrace::assertOrTerminate(const bool b,
                                                const ErrorReport e,
                                                const std::source_location &l) {
    if (!b) {
      this->terminate(e, l);
    }
  }  // end of assertOrTerminate

  template <typename ReturnType>
  ReturnType ErrorBacktrace::registerErrorMessage(
      const char *const e,
      const std::source_location &l)  //
      requires(internal::InvalidValueTraits<ReturnType>::isSpecialized) {
    std::ignore = this->registerErrorMessage(e, l);
    return internal::InvalidValueTraits<ReturnType>::getValue();
  }  // end of registerErrorMessage

  template <typename ReturnType>
  ReturnType ErrorBacktrace::registerErrorMessage(
      const ErrorReport e,
      const std::source_location &l) noexcept  //
      requires(internal::InvalidValueTraits<ReturnType>::isSpecialized) {
    std::ignore = this->registerErrorMessage(e, l);
    return internal::InvalidValueTraits<ReturnType>::getValue();
  }  // end of registerErrorMessage

#else

  inline void ErrorBacktrace::assertOrTerminate(const bool b, const char *msg) {
    if (!b) {
      this->terminate(msg);
    }
  }  // end of assertOrTerminate

  inline void ErrorBacktrace::assertOrTerminate(const bool b,
                                                const ErrorReport e) {
    if (!b) {
      this->terminate(e);
    }
  }  // end of assertOrTerminate

  template <typename ReturnType>
  ReturnType ErrorBacktrace::registerErrorMessage(const char *const e)  //
      requires(internal::InvalidValueTraits<ReturnType>::isSpecialized) {
    std::ignore = this->registerErrorMessage(e);
    return internal::InvalidValueTraits<ReturnType>::getValue();
  }  // end of registerErrorMessage

  template <typename ReturnType>
  ReturnType ErrorBacktrace::registerErrorMessage(
      const ErrorReport e) noexcept  //
      requires(internal::InvalidValueTraits<ReturnType>::isSpecialized) {
    std::ignore = this->registerErrorMessage(e);
    return internal::InvalidValueTraits<ReturnType>::getValue();
  }  // end of registerErrorMessage

#endif

}  // end of namespace mgis

#endif /* LIB_MGIS_ERRORBACKTRACE_IXX */
