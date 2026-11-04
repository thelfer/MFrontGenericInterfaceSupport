/*!
 * \file   MGIS/Context.ixx
 * \brief
 * \author Thomas Helfer
 * \date   24/02/2026
 */

#ifndef LIB_MGIS_CONTEXT_IXX
#define LIB_MGIS_CONTEXT_IXX

namespace mgis {

  template <typename... Args>
  std::ostream &Context::log(const VerbosityLevel l, Args &&...args) MGIS_NOEXCEPT {
    auto &os = this->log();
    if (this->getVerbosityLevel() >= l) {
      (os << ... << std::forward<Args>(args));
      os << std::endl;
    }
    return os;
  }  // end of log

  template <typename... Args>
  void Context::warning(Args &&...args) MGIS_NOEXCEPT {
    mgis::warning(this->log(), std::forward<Args>(args)...);
  }  // end of debug

  template <typename... Args>
  void Context::debug(Args &&...args) MGIS_NOEXCEPT {
    if (getVerbosityLevel() >= verboseDebug) {
      mgis::debug(this->log(), std::forward<Args>(args)...);
    }
  }  // end of debug

  inline Context::FailureHandler<Context::FailureHandlerPolicy::RAISE>
  Context::getThrowingFailureHandler() MGIS_NOEXCEPT {
    return this->getFailureHandler<FailureHandlerPolicy::RAISE>();
  }

  inline Context::FailureHandler<Context::FailureHandlerPolicy::ABORT>
  Context::getFatalFailureHandler() MGIS_NOEXCEPT {
    return this->getFailureHandler<FailureHandlerPolicy::ABORT>();
  }

}  // end of namespace mgis

#endif /* LIB_MGIS_CONTEXT_IXX */
