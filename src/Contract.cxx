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

  [[noreturn]] InvalidResult ContractChecker::registerErrorMessage(
      const char *const msg, const std::source_location &) {
    if constexpr (config::contract_violation_policy ==
                  config::ContractViolationPolicy::RAISE) {
      raise(msg);
    } else {
      ContractChecker::abort(msg);
    }
  }  // end of registerErrorMessage

#else

  [[noreturn]] InvalidResult ContractChecker::registerErrorMessage(
      const char* const msg) {
    if constexpr (config::contract_violation_policy ==
                  config::ContractViolationPolicy::RAISE) {
      raise(msg);
    } else {
      ContractChecker::abort(msg);
    }
  }  // end of registerErrorMessage

#endif

  void ContractChecker::abort(const char *const msg) noexcept {
    std::cerr << msg << '\n';
    std::abort();
  }  // end of abort

}  // namespace mgis
