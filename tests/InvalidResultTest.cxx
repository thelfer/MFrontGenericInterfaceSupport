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
  auto expect_false = [expect_true](const bool b) { expect_true(!b); };
  const auto r1 = []() -> std::optional<bool> { return InvalidResult{}; }();
  expect_false(r1.has_value());
  const auto r2 = []() -> std::optional<int> { return InvalidResult{}; }();
  expect_false(r2.has_value());
  const auto r3 = []() -> std::optional<size_type> {
    return InvalidResult{};
  }();
  expect_false(r3.has_value());
  const auto r4 = []() -> std::optional<real> { return InvalidResult{}; }();
  expect_false(r4.has_value());
  const auto r5 = []() -> std::optional<std::string> {
    return InvalidResult{};
  }();
  expect_false(r5.has_value());
  const auto r6 = []() -> std::unique_ptr<std::string> {
    return InvalidResult{};
  }();
  expect_true(r6.get() == nullptr);
  const auto r7 = []() -> std::shared_ptr<std::string> {
    return InvalidResult{};
  }();
  expect_true(r7.get() == nullptr);
  return EXIT_SUCCESS;
}
