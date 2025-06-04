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

  /*!
   * \brief a class used to report a contract violation
   *
   * The behaviour of this class is determined by the
   * `contract_violation_policy` variable.
   *
   * This class can be used in a `constexpr` context,
   * if no contrat violation is detected. If a
   * contract violation is detected, a compile-time
   * error is generated since `registerErrorMessage`
   * is not `constexpr`.
   */
  struct MGIS_EXPORT ContractViolationHandler final : AbstractErrorHandler {
    //! \brief destructor
    constexpr ContractViolationHandler() = default;
#ifdef MGIS_USE_SOURCE_LOCATION_INFORMATION
    /*!
     * \brief register a new error message
     * \param[in] e: error messgage
     */
    [[noreturn]] InvalidResult registerErrorMessage(
        const char *const,
        const std::source_location & =
            std::source_location::current()) override;
#else
    /*!
     * \brief register a new error message
     * \param[in] e: error message
     */
    [[noreturn]] InvalidResult registerErrorMessage(const char *const) override;
#endif
    //! \brief destructor
    constexpr ~ContractViolationHandler() override;

   protected:
    [[noreturn]] static void abort(const char *const) noexcept;
  };

  template <typename Type, typename... Args>
  inline constexpr auto is_check_preconditions_callable =
      requires(ContractViolationHandler e, Args... rargs) {
    { Type::checkPreconditions(e, rargs...) } -> std::same_as<bool>;
  };

  template <typename Type, typename... Args>
  constexpr void check_preconditions(Args &&...args) requires(
      is_check_preconditions_callable<Type, Args...>) {
    auto c = ContractViolationHandler{};
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
