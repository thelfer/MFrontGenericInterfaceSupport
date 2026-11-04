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
#include <vector>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MGIS/Utilities/Construct.hxx"

struct BuidableWithContext {
  BuidableWithContext(::mgis::Context &){};
};  // end of BuidableWithContext

struct BuidableWithAndWithoutContext {
  BuidableWithAndWithoutContext(::mgis::Context &) noexcept
      : build_with_context(true){};
  BuidableWithAndWithoutContext() noexcept {};
  bool build_with_context = false;
};  // end of BuidableWithContext

struct BuildableWithInitializerList {
  BuildableWithInitializerList(const std::size_t n, const double v)
      : values(n, v) {}
  BuildableWithInitializerList(std::initializer_list<double> l) : values(l) {}
  std::vector<double> values;
};  // end of BuildableWithInitializerList

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
    this->test6();
    this->test7();
    this->test8();
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

  void test6() {
    using namespace mgis;
    auto ctx = Context{};
    auto o = construct<BuidableWithContext>(ctx);
    TFEL_TESTS_ASSERT(isValid(o));
    auto po = make_unique<BuidableWithContext>(ctx);
    TFEL_TESTS_ASSERT(isValid(po));
  }  // end of test6

  void test7() {
    using namespace mgis;
    auto ctx = Context{};
    auto o = construct<BuidableWithAndWithoutContext>(ctx);
    TFEL_TESTS_ASSERT(isValid(o));
    TFEL_TESTS_CHECK_EQUAL(o->build_with_context, true);
    auto po = make_unique<BuidableWithAndWithoutContext>(ctx);
    TFEL_TESTS_ASSERT(isValid(po));
    TFEL_TESTS_CHECK_EQUAL(po->build_with_context, true);
  }

  void test8() {
    using namespace mgis;
    using Type = BuildableWithInitializerList;
    auto ctx = Context{};
    // the same constructor as for make_unique and make_shared must be called
    auto p1 = make_unique<Type>(ctx, std::size_t{3}, 1.5);
    auto p2 = make_unique_as<Type, Type>(ctx, std::size_t{3}, 1.5);
    auto p3 = make_shared<Type>(ctx, std::size_t{3}, 1.5);
    auto p4 = make_shared_as<Type, Type>(ctx, std::size_t{3}, 1.5);
    TFEL_TESTS_ASSERT((p1.get() != nullptr) && (p1->values.size() == 3));
    TFEL_TESTS_ASSERT((p2.get() != nullptr) && (p2->values.size() == 3));
    TFEL_TESTS_ASSERT((p3.get() != nullptr) && (p3->values.size() == 3));
    TFEL_TESTS_ASSERT((p4.get() != nullptr) && (p4->values.size() == 3));
  }  // end of test8
};

TFEL_TESTS_GENERATE_PROXY(ConstructTest, "ConstructTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main() {
  auto &m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("ConstructTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
