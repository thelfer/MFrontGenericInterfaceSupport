/*!
 * \file   MaterialPropertyTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   10/10/2022
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string_view>
#include "MGIS/Raise.hxx"
#include "MGIS/MaterialProperty/MaterialProperty.hxx"

static bool success = true;

static bool check(const bool b, const std::string_view msg) {
  if (!b) {
    success = false;
    std::cerr << msg << '\n';
  }
  return b;
}

static void test1(const std::string& library) {
  using namespace mgis;
  using namespace mgis::material_property;
  constexpr auto eps = real{1e-12};
  const auto mp = load(library, "Inconel600_YoungModulus");
  const auto op = OutOfBoundsPolicy::MGIS_MATERIALPROPERTY_NONE_POLICY;
  if (mp.inputs.size() != 1) {
    mgis::raise("invalid number of variables");
  }
  check(mp.output == "YoungModulus", "invalid output");
  check(mp.inputs[0] == "Temperature", "invalid first input");
  {
    constexpr auto T = real{500.93};
    constexpr auto Eref = real{203763998248.34};
    auto s = OutputStatus{};
    s.status = 2;
    const auto E = mp.fct(&s, &T, 1, op);
    check(s.status == 0,
          "invalid output status ('" + std::to_string(s.status) + "')");
    check(std::abs(E - Eref) < Eref * eps, "invalid output value");
  }
  {
    constexpr auto T = real{-1};
    auto s = OutputStatus{};
    s.status = 2;
    mp.fct(&s, &T, 1, op);
    check(s.status == -1,
          "invalid output status ('" + std::to_string(s.status) + "')");
    check(s.bounds_status == -1, "invalid bound status");
  }
}

static void test2(const std::string& library) {
  using namespace mgis;
  using namespace mgis::material_property;
  const auto mp = load(library, "MaterialPropertyBoundCheck");
  const auto op = OutOfBoundsPolicy::MGIS_MATERIALPROPERTY_NONE_POLICY;
  if (mp.inputs.size() != 1) {
    mgis::raise("invalid number of variables");
  }
  check(mp.output == "YoungModulus", "invalid output");
  check(mp.inputs[0] == "Temperature", "invalid first input");
  {
    constexpr auto T = real{250};
    auto s = OutputStatus{};
    s.status = 2;
    mp.fct(&s, &T, 1, op);
    check(s.status == 0,
          "invalid output status ('" + std::to_string(s.status) + "')");
  }
  {
    constexpr auto T = real{150};
    auto s = OutputStatus{};
    s.status = 2;
    mp.fct(&s, &T, 1, op);
    check(s.status == 0,
          "invalid output status ('" + std::to_string(s.status) + "')");
  }
  {
    constexpr auto T = real{500};
    auto s = OutputStatus{};
    s.status = 2;
    mp.fct(&s, &T, 1, op);
    check(s.status == 0,
          "invalid output status ('" + std::to_string(s.status) + "')");
  }
}

static void test3(const std::string& library) {
  using namespace mgis;
  using namespace mgis::material_property;
  constexpr auto eps = real{1e-12};
  const auto mp = load(library, "MaterialPropertyBoundCheck");
  const auto op = OutOfBoundsPolicy::MGIS_MATERIALPROPERTY_WARNING_POLICY;
  if (mp.inputs.size() != 1) {
    mgis::raise("invalid number of variables");
  }
  check(mp.output == "YoungModulus", "invalid output");
  check(mp.inputs[0] == "Temperature", "invalid first input");
  {
    constexpr auto T = real{250};
    auto s = OutputStatus{};
    s.status = 2;
    const auto E = mp.fct(&s, &T, 1, op);
    check(s.status == 0,
          "invalid output status ('" + std::to_string(s.status) + "')");
    check(std::abs(E - T) < T * eps, "invalid output value");
  }
  {
    constexpr auto T = real{150};
    auto s = OutputStatus{};
    s.status = 2;
    const auto E = mp.fct(&s, &T, 1, op);
    check(s.status == 1,
          "invalid output status ('" + std::to_string(s.status) + "')");
    check(s.bounds_status == 1, "invalid bound status");
    check(std::abs(E - T) < T * eps, "invalid output value");
  }
  {
    constexpr auto T = real{500};
    auto s = OutputStatus{};
    s.status = 2;
    const auto E = mp.fct(&s, &T, 1, op);
    check(s.status == 1, "invalid output status");
    check(s.status == 1,
          "invalid output status ('" + std::to_string(s.status) + "')");
    check(std::abs(E - T) < T * eps, "invalid output value");
  }
}

static void test4(const std::string& library) {
  using namespace mgis;
  using namespace mgis::material_property;
  constexpr auto eps = real{1e-12};
  const auto mp = load(library, "MaterialPropertyBoundCheck2");
  const auto op = OutOfBoundsPolicy::MGIS_MATERIALPROPERTY_WARNING_POLICY;
  if (mp.inputs.size() != 1) {
    mgis::raise("invalid number of variables");
  }
  check(mp.output == "YoungModulus", "invalid output");
  check(mp.inputs[0] == "Temperature", "invalid first input");
  {
    constexpr auto T = real{250};
    auto s = OutputStatus{};
    s.status = 2;
    const auto E = mp.fct(&s, &T, 1, op);
    check(s.status == 0,
          "invalid output status ('" + std::to_string(s.status) + "')");
    check(std::abs(E - T) < T * eps, "invalid output value");
  }
  {
    constexpr auto T = real{150};
    auto s = OutputStatus{};
    s.status = 2;
    const auto E = mp.fct(&s, &T, 1, op);
    check(s.status == 1,
          "invalid output status ('" + std::to_string(s.status) + "')");
    check(s.bounds_status == 1, "invalid bound status");
    check(std::abs(E - T) < T * eps, "invalid output value");
  }
  {
    constexpr auto T = real{500};
    auto s = OutputStatus{};
    s.status = 2;
    const auto E = mp.fct(&s, &T, 1, op);
    check(s.status == 0,
          "invalid output status ('" + std::to_string(s.status) + "')");
    check(std::abs(E - T) < T * eps, "invalid output value");
  }
}

static void test5(const std::string& library) {
  using namespace mgis;
  using namespace mgis::material_property;
  constexpr auto eps = real{1e-12};
  const auto mp = load(library, "MaterialPropertyBoundCheck3");
  const auto op = OutOfBoundsPolicy::MGIS_MATERIALPROPERTY_WARNING_POLICY;
  if (mp.inputs.size() != 1) {
    mgis::raise("invalid number of variables");
  }
  check(mp.output == "YoungModulus", "invalid output");
  check(mp.inputs[0] == "Temperature", "invalid first input");
  {
    constexpr auto T = real{250};
    auto s = OutputStatus{};
    s.status = 2;
    const auto E = mp.fct(&s, &T, 1, op);
    check(s.status == 0,
          "invalid output status ('" + std::to_string(s.status) + "')");
    check(std::abs(E - T) < T * eps, "invalid output value");
  }
  {
    constexpr auto T = real{150};
    auto s = OutputStatus{};
    s.status = 2;
    const auto E = mp.fct(&s, &T, 1, op);
    check(s.status == 0,
          "invalid output status ('" + std::to_string(s.status) + "')");
    check(std::abs(E - T) < T * eps, "invalid output value");
  }
  {
    constexpr auto T = real{500};
    auto s = OutputStatus{};
    s.status = 2;
    const auto E = mp.fct(&s, &T, 1, op);
    check(s.status == 1,
          "invalid output status ('" + std::to_string(s.status) + "')");
    check(s.bounds_status == 1, "invalid bound status");
    check(std::abs(E - T) < T * eps, "invalid output value");
  }
}

static void test6(const std::string& library) {
  using namespace mgis;
  using namespace mgis::material_property;
  const auto mp = load(library, "MaterialPropertyBoundCheck");
  const auto op = OutOfBoundsPolicy::MGIS_MATERIALPROPERTY_STRICT_POLICY;
  if (mp.inputs.size() != 1) {
    mgis::raise("invalid number of variables");
  }
  check(mp.output == "YoungModulus", "invalid output");
  check(mp.inputs[0] == "Temperature", "invalid first input");
  {
    constexpr auto T = real{250};
    auto s = OutputStatus{};
    s.status = 2;
    mp.fct(&s, &T, 1, op);
    check(s.status == 0,
          "invalid output status ('" + std::to_string(s.status) + "')");
  }
  {
    constexpr auto T = real{150};
    auto s = OutputStatus{};
    s.status = 2;
    mp.fct(&s, &T, 1, op);
    check(s.status == -1,
          "invalid output status ('" + std::to_string(s.status) + "')");
    check(s.bounds_status == -1, "invalid bound status");
  }
  {
    constexpr auto T = real{500};
    auto s = OutputStatus{};
    s.status = 2;
    mp.fct(&s, &T, 1, op);
    check(s.status == -1, "invalid output status");
    check(s.bounds_status == -1, "invalid bound status");
  }
}

static void test7(const std::string& library) {
  using namespace mgis;
  using namespace mgis::material_property;
  constexpr auto eps = real{1e-12};
  const auto mp = load(library, "MaterialPropertyErrnoCheck");
  const auto op = OutOfBoundsPolicy::MGIS_MATERIALPROPERTY_NONE_POLICY;
  if (mp.inputs.size() != 1) {
    mgis::raise("invalid number of variables");
  }
  check(mp.output == "YoungModulus", "invalid output");
  check(mp.inputs[0] == "Temperature", "invalid first input");
  {
    constexpr auto T = real{0.34};
    auto s = OutputStatus{};
    s.status = 2;
    const auto E = mp.fct(&s, &T, 1, op);
    check(s.status == 0,
          "invalid output status ('" + std::to_string(s.status) + "')");
    check(std::abs(E - std::acos(T)) < eps, "invalid output value");
  }
#ifndef __INTEL_LLVM_COMPILER
  {
    constexpr auto T = real{-2};
    auto s = OutputStatus{};
    s.status = 2;
    mp.fct(&s, &T, 1, op);
    check(s.status == -3,
          "invalid output status ('" + std::to_string(s.status) + "')");
  }
#endif /* __INTEL_LLVM_COMPILER */
}

static void test8(const std::string& library) {
  using namespace mgis;
  using namespace mgis::material_property;
  constexpr auto eps = real{1e-12};
  const auto mp = load(library, "MaterialPropertyErrnoCheck2");
  const auto op = OutOfBoundsPolicy::MGIS_MATERIALPROPERTY_NONE_POLICY;
  if (mp.inputs.size() != 1) {
    mgis::raise("invalid number of variables");
  }
  check(mp.output == "YoungModulus", "invalid output");
  check(mp.inputs[0] == "Temperature", "invalid first input");
  {
    constexpr auto T = real{0.34};
    auto s = OutputStatus{};
    s.status = 2;
    const auto E = mp.fct(&s, &T, 1, op);
    check(s.status == 0, "invalid output status");
    check(std::abs(E - std::acos(T)) < eps, "invalid output value");
  }
  {
    constexpr auto T = real{-2};
    auto s = OutputStatus{};
    s.status = 2;
    mp.fct(&s, &T, 1, op);
    check(s.status == -4, "invalid output status");
  }
}

int main(const int argc, const char* const* argv) {
  if (argc != 2) {
    std::cerr << "MaterialPropertyTest: invalid number of arguments\n";
    std::exit(-1);
  }
  try {
    test1(argv[1]);
    test2(argv[1]);
    test3(argv[1]);
    test4(argv[1]);
    test5(argv[1]);
    test6(argv[1]);
    test7(argv[1]);
    test8(argv[1]);
  } catch (std::exception& e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
