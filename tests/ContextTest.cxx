#include <vector>
#include <utility>
#include <iostream>
#include <string_view>
#include "MGIS/Context.hxx"

std::optional<std::vector<int>> f() { return std::vector<int>{1, 3, 6}; }

int main() {
  using namespace mgis;
  auto ctx = Context{};
  auto or_raise = ctx.getFailureHandler<>();
  auto v = or_raise(f());
  std::cout << v[0] << " " << v[1] << " " << v[2] << '\n';
  return EXIT_SUCCESS;
}
