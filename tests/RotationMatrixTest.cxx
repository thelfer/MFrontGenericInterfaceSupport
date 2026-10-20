/*!
 * \file   tests/RotationMatrixTest.cxx
 * \brief  tests of the size checks and of the copies of rotation matrices
 * \date   20/10/2026
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <vector>
#include <cstdlib>
#include "MGIS/Behaviour/RotationMatrix.hxx"

int main() {
  using namespace mgis;
  using namespace mgis::behaviour;
  auto expect_true = [](const bool b) {
    if (!b) {
      std::exit(-1);
    }
  };
  auto throws = [](auto f) {
    try {
      f();
    } catch (...) {
      return true;
    }
    return false;
  };
  const auto e1 = std::vector<real>{1, 0, 0, 1};
  const auto e2 = std::vector<real>{0, 1, 0, 0, 1, 0};
  const auto e3 = std::vector<real>{1, 0, 0, 0};
  // size checks
  expect_true(
      !throws([&] { RotationMatrix2D(e1, StorageMode::EXTERNAL_STORAGE); }));
  expect_true(throws([&] {
    RotationMatrix2D(std::span(e1).first(3), StorageMode::EXTERNAL_STORAGE);
  }));
  expect_true(!throws(
      [&] { RotationMatrix3D(e2, e2, StorageMode::EXTERNAL_STORAGE); }));
  expect_true(
      throws([&] { RotationMatrix3D(e3, e3, StorageMode::EXTERNAL_STORAGE); }));
  // a copy of a rotation matrix storing its values refers to its own values
  const auto r = RotationMatrix3D(e2, e2, StorageMode::LOCAL_STORAGE);
  const auto c = r;
  expect_true(c.a1.a.data() != r.a1.a.data());
  expect_true(c.a2.a.data() != r.a2.a.data());
  expect_true((c.a1.a.size() == 6u) && (c.a1.a[1] == 1));
  return EXIT_SUCCESS;
}
