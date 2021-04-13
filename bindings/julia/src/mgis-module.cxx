/*!
 * \file   bindings/julia/src/mgis-module.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/05/2019
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <jlcxx/jlcxx.hpp>

// #include <array>
// #include "MGIS/Julia/ArrayView.hxx"

// forward declarations
void declareThreadPool(jlcxx::Module& m);

// std::array<double, 3>& test() {
//   static std::array<double, 3> a = {1, 2, 3};
//   return a;
// }

JLCXX_MODULE define_mgis_module(jlcxx::Module& m) {
  //   m.method("test", []() {
  //     auto& a = test();
  //     return mgis::julia::ArrayView<double, 1>(a.data(), a.size());
  //   });
  declareThreadPool(m);
}
