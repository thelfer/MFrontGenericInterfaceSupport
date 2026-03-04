/*!
 * \file   tests/ConstructTest.cxx
 * \brief  This file contains unit test of the `construct`, `make_unique` and
 * `make_shared` functions
 * \author Thomas Helfer
 * \date   20/02/2026
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MGIS/Utilities/Construct.hxx"

struct ConstructTest final : public tfel::tests::TestCase {
  ConstructTest()
      : tfel::tests::TestCase("MGIS/Utilities", "ConstructTests") {
  }  // end of ConstructTest
  tfel::tests::TestResult execute() override {
    this->test1();
    this->test2();
    this->test3();
    this->test4();
    this->test5();
    return this->result;
  }  // end of execute

 private:
  struct Dummy {
    Dummy(const bool b) {
      if (b) {
        mgis::raise("invalid argument");
      }
    }
  };

  void test1() {
    using namespace mgis;
    auto ctx = Context{};
    auto d1 = construct<Dummy>(ctx, false);
    TFEL_TESTS_ASSERT(d1.has_value());
    TFEL_TESTS_ASSERT(ctx.empty());
    auto d2 = construct<Dummy>(ctx, true);
    TFEL_TESTS_ASSERT(!d2.has_value());
    TFEL_TESTS_ASSERT(!ctx.empty());
    TFEL_TESTS_CHECK_EQUAL(ctx.getRawErrorMessage(), "invalid argument");
    ctx.clearErrorMessages();
    auto d3 = MGIS_CONSTRUCT(Dummy, ctx, false);
    TFEL_TESTS_ASSERT(d3.has_value());
    TFEL_TESTS_ASSERT(ctx.empty());
    auto d4 = MGIS_CONSTRUCT(Dummy, ctx, true);
    TFEL_TESTS_ASSERT(!d4.has_value());
    TFEL_TESTS_ASSERT(!ctx.empty());
    TFEL_TESTS_CHECK_EQUAL(ctx.getRawErrorMessage(), "invalid argument");
  }  // end of test1

  void test2() {
    using namespace mgis;
    auto ctx = Context{};
    auto d1 = make_unique<Dummy>(ctx, false);
    TFEL_TESTS_ASSERT(d1.get() != nullptr);
    TFEL_TESTS_ASSERT(ctx.empty());
    auto d2 = make_unique<Dummy>(ctx, true);
    TFEL_TESTS_ASSERT(d2.get() == nullptr);
    TFEL_TESTS_ASSERT(!ctx.empty());
    TFEL_TESTS_CHECK_EQUAL(ctx.getRawErrorMessage(), "invalid argument");
    ctx.clearErrorMessages();
    auto d3 = MGIS_MAKE_UNIQUE(Dummy, ctx, false);
    TFEL_TESTS_ASSERT(d3.get() != nullptr);
    TFEL_TESTS_ASSERT(ctx.empty());
    auto d4 = MGIS_MAKE_UNIQUE(Dummy, ctx, true);
    TFEL_TESTS_ASSERT(d4.get() == nullptr);
    TFEL_TESTS_ASSERT(!ctx.empty());
    TFEL_TESTS_CHECK_EQUAL(ctx.getRawErrorMessage(), "invalid argument");
  }  // end of test2

  void test3() {
    using namespace mgis;
    auto ctx = Context{};
    auto d1 = make_shared<Dummy>(ctx, false);
    TFEL_TESTS_ASSERT(d1.get() != nullptr);
    TFEL_TESTS_ASSERT(ctx.empty());
    auto d2 = make_shared<Dummy>(ctx, true);
    TFEL_TESTS_ASSERT(d2.get() == nullptr);
    TFEL_TESTS_ASSERT(!ctx.empty());
    TFEL_TESTS_CHECK_EQUAL(ctx.getRawErrorMessage(), "invalid argument");
    ctx.clearErrorMessages();
    auto d3 = MGIS_MAKE_SHARED(Dummy, ctx, false);
    TFEL_TESTS_ASSERT(d3.get() != nullptr);
    TFEL_TESTS_ASSERT(ctx.empty());
    auto d4 = MGIS_MAKE_SHARED(Dummy, ctx, true);
    TFEL_TESTS_ASSERT(d4.get() == nullptr);
    TFEL_TESTS_ASSERT(!ctx.empty());
    TFEL_TESTS_CHECK_EQUAL(ctx.getRawErrorMessage(), "invalid argument");
  }  // end of test3

  void test4() {
    using namespace mgis;
    auto try_construct = [](const bool b) -> bool {
      auto ctx = Context{};
      MGIS_TRY_CONSTRUCT(Dummy, v, ctx, b);
      static_assert(std::is_same_v<decltype(v), Dummy &>,
                    "expected a non-const reference");
      std::ignore = v;
      return true;
    };
    auto try_make_unique = [](const bool b) -> bool {
      auto ctx = Context{};
      MGIS_TRY_MAKE_UNIQUE(Dummy, v, ctx, b);
      static_assert(std::is_same_v<decltype(*v), Dummy &>,
                    "expected a non-const reference");
      std::ignore = v;
      return true;
    };
    auto try_make_unique_as = [](const bool b) -> bool {
      auto ctx = Context{};
      MGIS_TRY_MAKE_UNIQUE_AS(Dummy, Dummy, v, ctx, b);
      static_assert(std::is_same_v<decltype(*v), Dummy &>,
                    "expected a non-const reference");
      std::ignore = v;
      return true;
    };
    auto try_make_shared = [](const bool b) -> bool {
      auto ctx = Context{};
      MGIS_TRY_MAKE_SHARED(Dummy, v, ctx, b);
      static_assert(std::is_same_v<decltype(*v), Dummy &>,
                    "expected a non-const reference");
      std::ignore = v;
      return true;
    };
    auto try_make_shared_as = [](const bool b) -> bool {
      auto ctx = Context{};
      MGIS_TRY_MAKE_SHARED_AS(Dummy, Dummy, v, ctx, b);
      static_assert(std::is_same_v<decltype(*v), Dummy &>,
                    "expected a non-const reference");
      std::ignore = v;
      return true;
    };
    TFEL_TESTS_ASSERT(try_construct(false));
    TFEL_TESTS_ASSERT(!try_construct(true));
    TFEL_TESTS_ASSERT(try_make_unique(false));
    TFEL_TESTS_ASSERT(!try_make_unique(true));
    TFEL_TESTS_ASSERT(try_make_unique_as(false));
    TFEL_TESTS_ASSERT(!try_make_unique_as(true));
    TFEL_TESTS_ASSERT(try_make_shared(false));
    TFEL_TESTS_ASSERT(!try_make_shared(true));
    TFEL_TESTS_ASSERT(try_make_shared_as(false));
    TFEL_TESTS_ASSERT(!try_make_shared_as(true));
  }

  void test5() {
    using namespace mgis;
    auto try_construct = [](const bool b) -> bool {
      auto ctx = Context{};
      MGIS_TRY_CONSTRUCT(const Dummy, v, ctx, b);
      static_assert(std::is_same_v<decltype(v), const Dummy &>,
                    "expected a const reference");
      std::ignore = v;
      return true;
    };
    auto try_make_unique = [](const bool b) -> bool {
      auto ctx = Context{};
      MGIS_TRY_MAKE_UNIQUE(const Dummy, v, ctx, b);
      static_assert(std::is_same_v<decltype(*v), const Dummy &>,
                    "expected a const reference");
      std::ignore = v;
      return true;
    };
    auto try_make_shared = [](const bool b) -> bool {
      auto ctx = Context{};
      MGIS_TRY_MAKE_SHARED(const Dummy, v, ctx, b);
      static_assert(std::is_same_v<decltype(*v), const Dummy &>,
                    "expected a const reference");
      std::ignore = v;
      return true;
    };
    TFEL_TESTS_ASSERT(try_construct(false));
    TFEL_TESTS_ASSERT(!try_construct(true));
    TFEL_TESTS_ASSERT(try_make_unique(false));
    TFEL_TESTS_ASSERT(!try_make_unique(true));
    TFEL_TESTS_ASSERT(try_make_shared(false));
    TFEL_TESTS_ASSERT(!try_make_shared(true));
  }
};

TFEL_TESTS_GENERATE_PROXY(ConstructTest, "ConstructTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto &m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("ConstructTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
