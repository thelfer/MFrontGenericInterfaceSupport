/*!
 * \file   src/Context.cxx
 * \brief  This file implements the `Context` class
 * \date   08/02/2023
 */

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

  Context::~Context() noexcept = default;

}  // end of namespace mgis
