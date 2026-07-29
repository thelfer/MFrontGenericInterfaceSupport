/*!
 * \file   bindings/python/src/Function.cxx
 * \brief
 * \author Thomas Helfer
 * \date   18/11/2025
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "MGIS/Python/NumPySupport.hxx"
#include "MGIS/Function/Function.hxx"
#include "MGIS/Function/BasicLinearSpace.hxx"

mgis::function::FunctionView<mgis::function::BasicLinearSpace,
                             mgis::function::FunctionDataLayoutDescription{},
                             false>
mgis_convert_to_span(const pybind11::array_t<double> &o) {
  using namespace mgis::function;
  const auto i = o.request();
  if (i.ndim != 1) {
    const auto s = static_cast<mgis::size_type>(i.size);
    const auto values = std::span{static_cast<const double *>(i.ptr), s};
    return FunctionView<BasicLinearSpace, FunctionDataLayoutDescription{},
                        false>{BasicLinearSpace{s}, values, 1};
  }
  mgis::raise("convert_to_span: expected one dimensional array");
}  // end of mgis_convert_to_span

void declareFunction(pybind11::module_ &m) {
  // FunctionView<BasicLinearSpace, FunctionDataLayoutDescription{}, false
  using mgis::real;
  //
  using BasicFunction =
      mgis::function::Function<mgis::function::BasicLinearSpace>;
  //
  using BasicFunctionView = mgis::function::FunctionView<
      mgis::function::BasicLinearSpace,
      mgis::function::FunctionDataLayoutDescription{}, true>;

  pybind11::class_<BasicFunction>(m, "BasicFunction",
                                  pybind11::buffer_protocol())
      .def_buffer([](BasicFunction &f) -> pybind11::buffer_info {
        if (f.getNumberOfComponents() == 1) {
          return pybind11::buffer_info(
              f.data().data(), sizeof(real),
              pybind11::format_descriptor<real>::format(), 1,
              {getSpaceSize(f.getSpace())}, {sizeof(real)});
        }
        return pybind11::buffer_info(
            f.data().data(), sizeof(real),
            pybind11::format_descriptor<real>::format(), 2,
            {getSpaceSize(f.getSpace()), f.getNumberOfComponents()},
            {sizeof(real) * f.getNumberOfComponents(), sizeof(real)});
      });

  pybind11::class_<BasicFunctionView>(m, "BasicFunctionView",
                                      pybind11::buffer_protocol())
      .def_buffer([](BasicFunctionView &f) -> pybind11::buffer_info {
        if (f.getNumberOfComponents() == 1) {
          return pybind11::buffer_info(
              f.data().data(), sizeof(real),
              pybind11::format_descriptor<real>::format(), 1,
              {getSpaceSize(f.getSpace())}, {sizeof(real) * f.getDataStride()});
        }
        return pybind11::buffer_info(
            f.data().data(), sizeof(real),
            pybind11::format_descriptor<real>::format(), 2,
            {getSpaceSize(f.getSpace()), f.getNumberOfComponents()},
            {sizeof(real) * f.getDataStride(), sizeof(real)});
      });

}  // end of declareFunction
