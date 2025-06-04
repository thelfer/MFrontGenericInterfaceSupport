/*!
 * \file   MGIS/Contract.hxx
 * \brief
 * \author Thomas Helfer
 * \date   03/06/2025
 */

#ifndef LIB_MGIS_CONTRACT_HXX
#define LIB_MGIS_CONTRACT_HXX

#include "MGIS/Raise.hxx"
#include "MGIS/Config.hxx"
#include "MGIS/InvalidResult.hxx"
#include "MGIS/AbstractErrorHandler.hxx"

namespace mgis {

  struct MGIS_EXPORT ContractViolationReporter final : AbstractErrorHandler {
    constexpr ContractViolationReporter() = default;
#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
    [[noreturn]] InvalidResult registerErrorMessage(
        const char *const,
        const std::source_location & =
            std::source_location::current()) override;
#else
    [[noreturn]] InvalidResult registerErrorMessage(const char *const) override;
#endif
    constexpr ~ContractViolationReporter() override;

   protected:
    [[noreturn]] static void abort(const char *const) noexcept;
  };

  template <typename Type, typename... Args>
  inline constexpr auto is_check_preconditions_callable =
      requires(ContractViolationReporter e, Args... rargs) {
    { Type::checkPreconditions(e, rargs...) } -> std::same_as<bool>;
  };

  template <typename Type, typename... Args>
  constexpr void check_preconditions(Args &&...args) requires(
      is_check_preconditions_callable<Type, Args...>) {
    auto c = ContractViolationReporter{};
    static_cast<void>(Type::checkPreconditions(c, std::forward<Args>(args)...));
  }  // end of check_preconditions

  template <bool b>
  struct PreconditionsCheck {};

  inline constexpr auto no_precondition_check = PreconditionsCheck<false>{};
  inline constexpr auto preconditions_check = PreconditionsCheck<true>{};

  template <typename Child>
  struct PreconditionsChecker {
    PreconditionsChecker(PreconditionsChecker &&) = default;
    PreconditionsChecker(const PreconditionsChecker &) = default;
    template <bool b, typename... Args>
    constexpr PreconditionsChecker(const PreconditionsCheck<b> &,
                                   Args &&...args) {
      if constexpr (b) {
        check_preconditions<Child>(std::forward<Args>(args)...);
      }
    }
  };

}  // namespace mgis

#include "MGIS/Contract.ixx"

#endif /* LIB_MGIS_CONTRACT_HXX */
