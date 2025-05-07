/*!
 * \file   src/ErrorBacktrace.cxx
 * \brief  This file implements the `ErrorBacktrace` class.
 * \date   04/11/2022
 */

#include <sstream>
#include <iostream>
#include "MGIS/ErrorBacktrace.hxx"

namespace mgis::internal {

  //! \return a flag stating if registring an error is fatal, i.e. exits by
  //! calling `std::exit`.
  static bool &getErrorReportingPolicy() {
    static bool isErrorReportingFatal = false;
    return isErrorReportingFatal;
  }  // end of getErrorReportingPolicy

}  // end of namespace mgis::internal

namespace mgis {

  static std::vector<std::string> split(const std::string &s) {
    auto result = std::vector<std::string>{};
    auto ss = std::stringstream{s};
    for (std::string line; std::getline(ss, line, '\n');) {
      result.push_back(line);
    }
    return result;
  }  // end of split

  void ErrorBacktrace::setErrorReportingAsFatal() noexcept {
    ::mgis::internal::getErrorReportingPolicy() = true;
  }  // end of setIfErrorReportingIsFatal

  void ErrorBacktrace::unsetErrorReportingAsFatal() noexcept {
    ::mgis::internal::getErrorReportingPolicy() = false;
  }  // end of setIfErrorReportingIsFatal

#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
  InvalidResult ErrorBacktrace::registerErrorMessage(
      const ErrorReport e, const std::source_location &l) noexcept {
    try {
      this->errorMessages_.push_back(
          {l.file_name(), l.function_name(), l.line(), e});
    } catch (...) {
    }
    this->treatFatalCase_();
    return {};
  }  // end of registerErrorMessage
#else
  InvalidResult ErrorBacktrace::registerErrorMessage(
      const ErrorReport e) noexcept {
    return this->registerErrorMessageWithoutSourceLocation(e);
  }  // end of registerErrorMessage
#endif

  InvalidResult ErrorBacktrace::registerErrorMessageWithoutSourceLocation(
      const ErrorReport e) noexcept {
    try {
      auto em = ErrorMessage{};
      em.e = e;
      this->errorMessages_.push_back(em);
    } catch (...) {
    }
    this->treatFatalCase_();
    return {};
  }  // end of registerErrorMessage

  std::string ErrorBacktrace::getErrorMessage(const ErrorReport &e) noexcept {
    if (std::holds_alternative<const char *>(e)) {
      return std::get<const char *>(e);
    }
    if (std::holds_alternative<std::string>(e)) {
      return std::get<std::string>(e);
    }
    const auto &er = std::get<std::pair<int, ErrorReportFunction>>(e);
    return er.second(er.first);
  }  // end of ErrorBacktrace::getErrorMessage

  std::string ErrorBacktrace::getErrorMessage(const bool b) noexcept {
#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
    auto msg = this->getErrorMessage_();
    if (b) {
      this->clearErrorMessages();
    }
    return msg;
#else
    return this->getRawErrorMessage(b);
#endif
  }  // end of getErrorMessage

  std::string ErrorBacktrace::getRawErrorMessage(const bool b) noexcept {
    auto msg = this->getRawErrorMessage_();
    if (b) {
      this->clearErrorMessages();
    }
    return msg;
  }  // end of getRawErrorMessage

  std::string ErrorBacktrace::getErrorMessage_() const noexcept {
#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
    auto msg = std::string{};
    auto prefix = std::string{};
    try {
      auto p = this->errorMessages_.rbegin();
      const auto pe = this->errorMessages_.rend();
      for (; p != pe; ++p) {
        auto firstLine = true;
        for (const auto &l : split(ErrorBacktrace::getErrorMessage(p->e))) {
          if (!msg.empty()) {
            msg += '\n';
          }
          msg += prefix;
          if (firstLine) {
            if (p->file_name != nullptr) {
              msg += p->file_name;
              msg += ":" + std::to_string(p->line) + ": in function '";
              msg += p->function_name;
              msg += "': ";
            }
          }
          msg += l;
          firstLine = false;
        }
        if (prefix.empty()) {
          prefix = ' ';
        }
        prefix = '*' + prefix;
      }
    } catch (...) {
    }
    return msg;
#else
    return this->getRawErrorMessage_();
#endif
  }  // end of getErrorMessage_

  std::string ErrorBacktrace::getRawErrorMessage_() const noexcept {
    auto msg = std::string{};
    auto prefix = std::string{};
    try {
      auto p = this->errorMessages_.rbegin();
      const auto pe = this->errorMessages_.rend();
      for (; p != pe; ++p) {
        for (const auto &l : split(ErrorBacktrace::getErrorMessage(p->e))) {
          if (!msg.empty()) {
            msg += '\n';
          }
          msg += prefix;
          msg += l;
        }
        if (prefix.empty()) {
          prefix = ' ';
        }
        prefix = '*' + prefix;
      }
    } catch (...) {
    }
    return msg;
  }  // end of getRawErrorMessage_

  bool ErrorBacktrace::empty() const noexcept {
    return this->errorMessages_.empty();
  }  // end of empty

  void ErrorBacktrace::clearErrorMessages() noexcept {
    return this->errorMessages_.clear();
  }  // end of clearErrorMessages

  void ErrorBacktrace::treatFatalCase_() const noexcept {
    if (::mgis::internal::getErrorReportingPolicy()) {
      std::cerr << this->getErrorMessage_() << std::endl;
      std::exit(-1);
    }
  }  // end of treatFatalCase_

  ErrorBacktrace::~ErrorBacktrace() noexcept = default;

#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
  InvalidResult registerExceptionInErrorBacktrace(
      ErrorBacktrace &e, const std::source_location &l) noexcept {
    try {
      throw;
    } catch (std::exception &exception) {
      e.registerErrorMessage(std::string{exception.what()}, l);
    } catch (...) {
      e.registerErrorMessage("unknown exception thrown");
    }
    return {};
  }
#else
  InvalidResult registerExceptionInErrorBacktrace(ErrorBacktrace &e) noexcept {
    return registerExceptionInErrorBacktraceWithoutSourceLocation(e);
  }
#endif

  InvalidResult registerExceptionInErrorBacktraceWithoutSourceLocation(
      ErrorBacktrace &e) noexcept {
    try {
      throw;
    } catch (std::exception &exception) {
      e.registerErrorMessageWithoutSourceLocation(
          std::string{exception.what()});
    } catch (...) {
      e.registerErrorMessageWithoutSourceLocation("unknown exception thrown");
    }
    return {};
  }

}  // namespace mgis
