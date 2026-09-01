#include <string>
#include <cstdlib>
#include "MGIS/Config.hxx"
#include "MGIS/InvalidResult.hxx"

int main() {
  using namespace mgis;
  auto expect_true = [](const bool b) {
    if (!b) {
      std::exit(-1);
    }
  };
  //
  auto expect_false = [expect_true](const bool b) { expect_true(!b); };
  const auto r1 = []() -> std::optional<bool> { return InvalidResult{}; }();
  expect_false(r1.has_value());
  expect_true(isInvalid(r1));
  expect_false(isValid(r1));
  //
  const auto r2 = []() -> std::optional<int> { return InvalidResult{}; }();
  expect_false(r2.has_value());
  expect_true(isInvalid(r2));
  expect_false(isValid(r2));
  //
  const auto r3 = []() -> std::optional<size_type> {
    return InvalidResult{};
  }();
  expect_false(r3.has_value());
  expect_true(isInvalid(r3));
  expect_false(isValid(r3));
  //
  const auto r4 = []() -> std::optional<real> { return InvalidResult{}; }();
  expect_false(r4.has_value());
  expect_true(isInvalid(r4));
  expect_false(isValid(r4));
  //
  const auto r5 = []() -> std::optional<std::string> {
    return InvalidResult{};
  }();
  expect_false(r5.has_value());
  expect_true(isInvalid(r5));
  expect_false(isValid(r5));
  //
  const auto r6 = []() -> std::unique_ptr<std::string> {
    return InvalidResult{};
  }();
  expect_true(r6.get() == nullptr);
  expect_true(isInvalid(r6));
  expect_false(isValid(r6));
  //
  const auto r7 = []() -> std::shared_ptr<std::string> {
    return InvalidResult{};
  }();
  expect_true(r7.get() == nullptr);
  expect_true(isInvalid(r7));
  expect_false(isValid(r7));
  //
  InvalidResult r8;
  expect_true(isInvalid(r8));
  expect_false(isValid(r8));
  return EXIT_SUCCESS;
}
