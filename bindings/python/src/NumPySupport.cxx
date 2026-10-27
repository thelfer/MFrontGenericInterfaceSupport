/*!
 * \file   NumPySupport.cxx
 * \brief
 * \author Thomas Helfer
 * \date   08/11/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <iostream>
#include "MGIS/Raise.hxx"
#include "MGIS/Python/NumPySupport.hxx"

namespace mgis::python {

  pybind11::array_t<double> wrapInNumPyArray(std::span<double>& v) {
    return pybind11::array_t<double>({v.size()}, {sizeof(double)}, v.data(),
                                     pybind11::none());
  }  // end of wrapInNumPyArray

  pybind11::array_t<double> wrapInNumPyArray(std::vector<double>& v) {
    auto view = std::span<double>(v);
    return wrapInNumPyArray(view);
  }  // end of wrapInNumPyArray

  pybind11::array_t<double> wrapInNumPyArray(
      std::variant<std::span<double>, std::vector<double>>& v) {
    if (std::holds_alternative<std::vector<double>>(v)) {
      return wrapInNumPyArray(std::get<std::vector<double>>(v));
    }
    return wrapInNumPyArray(std::get<std::span<double>>(v));
  }  // end of wrapInNumPyArray

  pybind11::array_t<double> wrapInNumPyArray(std::span<double>& v,
                                             const size_type nc) {
    if (nc == 0) {
      return pybind11::array_t<double>();
    }
    return pybind11::array_t<double>(
        {v.size() / nc, nc}, /* Buffer dimensions */
        {nc * sizeof(double), sizeof(double)}, v.data(), pybind11::none());
  }  // end of wrapInNumPyArray

  pybind11::array_t<double> wrapInNumPyArray(std::vector<double>& v,
                                             const size_type nc) {
    auto view = std::span<double>(v);
    return wrapInNumPyArray(view, nc);
  }  // end of wrapInNumPyArray

  pybind11::array_t<double> wrapInNumPyArray(
      std::variant<std::span<double>, std::vector<double>>& v,
      const size_type nc) {
    if (std::holds_alternative<std::vector<double>>(v)) {
      return wrapInNumPyArray(std::get<std::vector<double>>(v), nc);
    }
    return wrapInNumPyArray(std::get<std::span<double>>(v), nc);
  }  // end of wrapInNumPyArray

  pybind11::array_t<double> wrapInNumPyArray(std::span<double>& v,
                                             const size_type nl,
                                             const size_type nc) {
    if ((nc == 0) || (nl == 0)) {
      return pybind11::array_t<double>();
    }
    return pybind11::array_t<double>(
        {v.size() / (nc * nl), nl, nc},
        {(nc * nl) * sizeof(double), nc * sizeof(double), sizeof(double)},
        v.data(), pybind11::none());
  }  // end of wrapInNumPyArray

  pybind11::array_t<double> wrapInNumPyArray(std::vector<double>& v,
                                             const size_type nl,
                                             const size_type nc) {
    auto view = std::span<double>(v);
    return wrapInNumPyArray(view, nl, nc);
  }  // end of wrapInNumPyArray

  pybind11::array_t<double> wrapInNumPyArray(
      std::variant<std::span<double>, std::vector<double>>& v,
      const size_type nl,
      const size_type nc) {
    if (std::holds_alternative<std::vector<double>>(v)) {
      return wrapInNumPyArray(std::get<std::vector<double>>(v), nl, nc);
    }
    return wrapInNumPyArray(std::get<std::span<double>>(v), nl, nc);
  }  // end of wrapInNumPyArray

  std::span<mgis::real> mgis_convert_to_span(
      const pybind11::array_t<double>& o) {
    const auto i = o.request();
    if (i.ndim != 1) {
      mgis::raise("convert_to_span: expected one dimensional array");
    }
    if ((i.size > 1) &&
        (i.strides[0] != static_cast<pybind11::ssize_t>(sizeof(double)))) {
      mgis::raise("convert_to_span: expected a contiguous array");
    }
    return {static_cast<double*>(i.ptr), static_cast<mgis_size_type>(i.size)};
  }  // end of mgis_convert_to_span

  std::span<mgis::real> mgis_convert_to_mutable_span(
      const pybind11::object& o) {
    const auto a = pybind11::array_t<double>(o);
    // the values written in a converted copy would be lost
    if ((!a.is(o)) && (!pybind11::module_::import("numpy")
                            .attr("may_share_memory")(a, o)
                            .cast<bool>())) {
      mgis::raise(
          "convert_to_span: expected an array of float64 values "
          "(the values would be written in a converted copy)");
    }
    if (!a.writeable()) {
      mgis::raise("convert_to_span: expected a writable array");
    }
    return mgis_convert_to_span(a);
  }  // end of mgis_convert_to_mutable_span

}  // end of namespace mgis::python
