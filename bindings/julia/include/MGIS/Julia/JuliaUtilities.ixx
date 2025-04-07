/*!
 * \file   JuliaUtilities.ixx
 * \brief
 * \author th202608
 * \date   16/05/2019
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_JULIA_JULIAUTILITIES_IXX
#define LIB_MGIS_JULIA_JULIAUTILITIES_IXX

#include <cstdint>
#include <algorithm>

namespace mgis {

  namespace julia {

    template <typename T>
    void expose_std_vector(jlcxx::Module& m, const char* const n) {
      m.add_type<std::vector<T>>(n)
          .method("length",
                  [](const std::vector<T>& v) noexcept -> std::int64_t {
                    return static_cast<std::int64_t>(v.size());
                  })
          .method("setindex!",
                  [](std::vector<T>& v, const T& value, const std::int64_t i) {
                    if (i == 0) {
                      mgis::raise<std::range_error>("invalid index");
                    }
                    v.at(i - 1) = value;
                  })
          .method("getindex", [](const std::vector<T>& v,
                                 const std::int64_t i) {
            if (i <= 0) {
              mgis::raise<std::range_error>("invalid index");
            }
            return v.at(static_cast<typename std::vector<T>::size_type>(i - 1));
          });
      //           .method("=", [](std::vector<T>& v, const jlcxx::ArrayRef<T>&
      //           s) {
      //             assign(v, s);
      //           });
    }  // end of expose_std_vector

    template <typename T>
    jlcxx::ConstArray<T, 1u> make_view(const std::vector<T>& v) {
      return jlcxx::make_const_array(v.data(), v.size());
    }  // end of make_view

    template <typename T>
    void assign(std::vector<T>& d, const jlcxx::ArrayRef<T>& s) {
      if (d.size() != s.size()) {
        mgis::raise<std::range_error>(
            "assign: both objects have unmatched size");
      }
      for (decltype(d.size()) i = 0; i != d.size(); ++i) {
        d[i] = s[i];
      }
    }  // end of assign

  }  // end of namespace julia

}  // end of namespace mgis

#endif /* LIB_MGIS_JULIA_JULIAUTILITIES_IXX */
