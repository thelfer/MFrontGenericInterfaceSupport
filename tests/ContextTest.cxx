#include <vector>
#include <utility>
#include <iostream>
#include <string_view>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MGIS/Context.hxx"

struct ContextTest final : public tfel::tests::TestCase {
  ContextTest()
      : tfel::tests::TestCase("MGIS", "ContextTests") {}  // end of ContextTest
  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
    this->test3();
    this->test4();
    return this->result;
  }

 private:
  [[nodiscard]] static std::optional<std::vector<int>> f() {
    return std::vector<int>{1, 3, 6};
  }
  [[nodiscard]] static std::shared_ptr<std::vector<int>> f2() {
    return std::make_shared<std::vector<int>>(std::vector<int>{1, 3, 6});
  }
  [[nodiscard]] static std::unique_ptr<std::vector<int>> f3() {
    return std::make_unique<std::vector<int>>(std::vector<int>{1, 3, 6});
  }
  template <typename T>
  static constexpr bool lvalue_assign = requires(mgis::Context& c,
                                                 const std::optional<T>& v3) {
    v3 | c.getFailureHandler<>();
  };
  void test1() {
    using namespace mgis;
    auto ctx = Context{};
    auto or_raise = ctx.getFailureHandler<>();
    auto v = or_raise(f());
    TFEL_TESTS_ASSERT(v.size() == 3);
    TFEL_TESTS_CHECK_EQUAL(v[0], 1);
    TFEL_TESTS_CHECK_EQUAL(v[1], 3);
    TFEL_TESTS_CHECK_EQUAL(v[2], 6);
    auto v2 = f() | or_raise;
    TFEL_TESTS_ASSERT(v2.size() == 3);
    TFEL_TESTS_CHECK_EQUAL(v2[0], 1);
    TFEL_TESTS_CHECK_EQUAL(v2[1], 3);
    TFEL_TESTS_CHECK_EQUAL(v2[2], 6);
    TFEL_TESTS_CHECK_THROW(std::optional<int>() | or_raise, std::runtime_error);
    //
    TFEL_TESTS_STATIC_ASSERT(!lvalue_assign<int>);
  }  // end of test1
  void test2() {
    using namespace mgis;
    auto ctx = Context{};
    auto or_raise = ctx.getFailureHandler<>();
    int v = 12;
    auto& v2 = OptionalReference<int>(&v) | or_raise;
    TFEL_TESTS_CHECK_EQUAL(&v2, &v);
  }  // end of test2
  void test3() {
    using namespace mgis;
    auto copy = [](const auto v) { return v; };
    auto ctx = Context{};
    auto or_raise = ctx.getFailureHandler<>();
    auto v = or_raise(f2());
    TFEL_TESTS_ASSERT(v->size() == 3);
    TFEL_TESTS_CHECK_EQUAL(v->at(0), 1);
    TFEL_TESTS_CHECK_EQUAL(v->at(1), 3);
    TFEL_TESTS_CHECK_EQUAL(v->at(2), 6);
    auto v2 = copy(v) | or_raise;
    TFEL_TESTS_CHECK_EQUAL(v.get(), v2.get());
    TFEL_TESTS_CHECK_THROW(std::shared_ptr<int>() | or_raise,
                           std::runtime_error);
  }
  void test4() {
    using namespace mgis;
    auto ctx = Context{};
    auto or_raise = ctx.getFailureHandler<>();
    auto v = or_raise(f3());
    TFEL_TESTS_ASSERT(v->size() == 3);
    TFEL_TESTS_CHECK_EQUAL(v->at(0), 1);
    TFEL_TESTS_CHECK_EQUAL(v->at(1), 3);
    TFEL_TESTS_CHECK_EQUAL(v->at(2), 6);
    TFEL_TESTS_CHECK_THROW(std::unique_ptr<int>() | or_raise,
                           std::runtime_error);
  }
};

TFEL_TESTS_GENERATE_PROXY(ContextTest, "ContextTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("ContextTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
