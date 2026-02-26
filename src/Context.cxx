/*!
 * \file   src/Context.cxx
 * \brief  This file implements the `Context` class
 * \date   08/02/2023
 */

#include <fstream>
#include "MGIS/LogStream.hxx"
#include "MGIS/Context.hxx"

namespace mgis {

  Context::Context() noexcept
      : verbosity(mgis::getDefaultVerbosityLevel()) {  //
  }                                                    // end of Context

  Context::Context(const ContextInitializer &i) noexcept
      : verbosity(i.verbosity) {}  // end of Context

  const VerbosityLevel &Context::getVerbosityLevel() const noexcept {
    return this->verbosity;
  }  // end of getVerbosityLevel

  void Context::setVerbosityLevel(const VerbosityLevel l) noexcept {
    this->verbosity = l;
  }  // end of setVerbosityLevel

  void Context::setLogStream(std::ostream &s) noexcept {
    this->log_stream = &s;
  }  // end of setLogStream

  void Context::setLogStream(std::shared_ptr<std::ostream> s) noexcept {
    if (s.get() == nullptr) {
      this->resetLogStream();
      return;
    }
    this->log_stream = s;
  }  // end of setLogStream

  void Context::resetLogStream() noexcept {
    this->log_stream = std::monostate{};
  }

  void Context::disableLogStream() noexcept {
    /*!
     * \brief a buffer which allows to create no-op output streams, i.e. streams
     * that does not nothing.
     *
     * see
     * https://stackoverflow.com/questions/11826554/standard-no-op-output-stream
     * for details
     */
    struct NullBuffer : public std::streambuf {
      int overflow(int c) override { return c; }
    };
    static NullBuffer null_buffer;
    static std::ostream null_stream(&null_buffer);
    this->setLogStream(null_stream);
  }

  std::ostream &Context::log() noexcept {
    if (std::holds_alternative<std::ostream *>(this->log_stream)) {
      auto *const ptr = std::get<std::ostream *>(this->log_stream);
      if (ptr == nullptr) {
        return ::mgis::getDefaultLogStream();
      }
      return *ptr;
    } else if (std::holds_alternative<std::shared_ptr<std::ostream>>(
                   this->log_stream)) {
      auto ptr = std::get<std::shared_ptr<std::ostream>>(this->log_stream);
      if (ptr.get() == nullptr) {
        return ::mgis::getDefaultLogStream();
      }
      return *ptr;
    }
    return ::mgis::getDefaultLogStream();
  }  // end of log

  void Context::abort() {
    const auto msg = this->getErrorMessage();
    if (!msg.empty()) {
      this->log() << msg << std::endl;
    } else {
      this->log() << "fatal error" << std::endl;
    }
    std::abort();
  }  // end of abort

  Context::~Context() noexcept = default;

}  // end of namespace mgis
