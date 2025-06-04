/*!
 * \file   Contract.cxx
 * \brief
 * \author Thomas Helfer
 * \date   03/06/2025
 */

#include <cstdlib>
#include <iostream>
#include "MGIS/Contract.hxx"

namespace mgis {

#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION

  [[noreturn]] InvalidResult ContractViolationHandler::registerErrorMessage(
      const char *const msg, const std::source_location &) {
    if constexpr (config::contract_violation_policy ==
                  config::ContractViolationPolicy::RAISE) {
      raise(msg);
    } else {
      ContractViolationHandler::abort(msg);
    }
  }  // end of registerErrorMessage

#else

  [[noreturn]] InvalidResult ContractViolationHandler::registerErrorMessage(
      const char* const msg) {
    if constexpr (config::contract_violation_policy ==
                  config::ContractViolationPolicy::RAISE) {
      raise(msg);
    } else {
      ContractViolationHandler::abort(msg);
    }
  }  // end of registerErrorMessage

#endif

  void ContractViolationHandler::abort(const char *const msg) noexcept {
    std::cerr << msg << '\n';
    std::abort();
  }  // end of abort

}  // namespace mgis
