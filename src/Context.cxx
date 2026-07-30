/*!
 * \file   src/Context.cxx
 * \brief  This file implements the `Context` class
 * \date   08/02/2023
 */

#include <fstream>
#include <algorithm>
#include "MGIS/LogStream.hxx"
#include "MGIS/Context.hxx"

namespace mgis {

  Context::Context() noexcept : verbosity(mgis::getDefaultVerbosityLevel()) {
    // Initialisation de la racine du profilage
    this->root_profiling_data.name = "root";
    this->root_profiling_data.calls = 1;
    this->profiling_stack.push_back(&this->root_profiling_data);
  }  // end of Context

  Context::Context(const ContextInitializer &i) noexcept
      : verbosity(i.verbosity) {
    // Initialisation de la racine du profilage
    this->root_profiling_data.name = "root";
    this->root_profiling_data.calls = 1;
    this->profiling_stack.push_back(&this->root_profiling_data);
  }  // end of Context

  const VerbosityLevel &Context::getVerbosityLevel() const noexcept {
    return this->verbosity;
  }  // end of getVerbosityLevel

  void Context::setVerbosityLevel(const VerbosityLevel l) noexcept {
    this->verbosity = l;
  }  // end of setVerbosityLevel

  void Context::enableProfiling(const bool b) noexcept {
    this->profiling_enabled = b;
  }  // end of enableProfiling

  bool Context::isProfilingEnabled() const noexcept {
    return this->profiling_enabled;
  }  // end of isProfilingEnabled

  ProfilingSection Context::startNewProfiling(std::string name,
                                              bool enabled) noexcept {
    return ProfilingSection{*this, std::move(name), enabled};
  }  // end of startNewProfiling

  void Context::pushProfilingNode(std::string name) noexcept {
    if (this->profiling_stack.empty()) return;

    ProfilingData *current = this->profiling_stack.back();

    auto it =
        std::find_if(current->children.begin(), current->children.end(),
                     [&name](const std::unique_ptr<ProfilingData> &child) {
                       return child->name == name;
                     });

    if (it != current->children.end()) {
      this->profiling_stack.push_back(it->get());
    } else {
      auto new_node = std::make_unique<ProfilingData>();
      new_node->name = std::move(name);

      ProfilingData *ptr = new_node.get();
      current->children.push_back(std::move(new_node));
      this->profiling_stack.push_back(ptr);
    }
  }  // end of pushProfilingNode

  void Context::popProfilingNode(double dt) noexcept {
    if (this->profiling_stack.size() > 1) {
      ProfilingData *current = this->profiling_stack.back();
      current->time_in_seconds += dt;
      current->calls += 1;

      this->profiling_stack.pop_back();
    } else if (this->profiling_stack.size() == 1) {
      this->profiling_stack.back()->time_in_seconds += dt;
    }
  }  // end of popProfilingNode

  const ProfilingData &Context::getProfilingResultTree() const noexcept {
    return this->root_profiling_data;
  }  // end of getProfilingResultTree

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

  std::shared_ptr<std::ostream> Context::getLogStreamPointer() const noexcept {
    if (std::holds_alternative<std::shared_ptr<std::ostream>>(
            this->log_stream)) {
      return std::get<std::shared_ptr<std::ostream>>(this->log_stream);
    }
    return {};
  }  // end of getLogStreamPointer

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