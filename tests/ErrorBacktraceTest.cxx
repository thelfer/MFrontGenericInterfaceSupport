#include <string_view>
#include "MGIS/Raise.hxx"
#include "MGIS/ErrorBacktrace.hxx"

namespace mgis::internal {

  void processNext(){};

  bool f3(ErrorBacktrace &e) { return e.registerErrorMessage("invalid call"); }

  bool f2(ErrorBacktrace &e) {
    if (!f3(e)) {
      // f2 fails, but we don't have any more information
      // to add, so we just return
      return false;
    }
    processNext();
    return true;
  }

  bool f1(ErrorBacktrace &e) {
    if (!f2(e)) {
      return e.registerErrorMessage("invalid call to f2");
    }
    return true;
  }

}  // end of namespace mgis::internal

int main() {
  using namespace mgis;
  auto expect_true = [](const bool b) {
    if (!b) {
      std::exit(-1);
    }
  };
  auto expect_false = [expect_true](const bool b) { expect_true(!b); };
  auto expect_eq = [expect_true](const std::string_view a,
                                 const std::string_view b) {
    expect_true(a == b);
  };
  {
    ErrorBacktrace e;
    try {
      raise("invalid argument");
    } catch (...) {
      registerExceptionInErrorBacktrace(e);
    }
    expect_eq(e.getRawErrorMessage(), "invalid argument");
  }  // end of ConstructTest
  {
    ErrorBacktrace e;
    expect_false(::mgis::internal::f1(e));
    expect_eq(e.getRawErrorMessage(), "invalid call to f2\n* invalid call");
  }  // end of ConstructTest
  {
    ErrorBacktrace e;
    e.registerErrorMessage(
        "first error message\nwith details on multiple lines");
    e.registerErrorMessage("parent error message");
    expect_eq(e.getRawErrorMessage(),
              "parent error message\n"
              "* first error message\n"
              "* with details on multiple lines");
  }  // end of ConstructTest
  return EXIT_SUCCESS;
}
