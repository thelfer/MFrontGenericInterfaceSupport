/*!
 * \file   bindings/python/src/Model.cxx
 * \brief
 * \author Thomas Helfer
 * \date   06/11/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "MGIS/Model/Model.hxx"

// forward declaration
void declareModel(pybind11::module_& m);

static mgis::model::Model Model_loadFromDatabase(
    const std::string& n,
    const mgis::behaviour::Hypothesis h,
    const std::optional<std::string>& m) {
  return mgis::model::loadFromDatabase(
      {.name = n, .hypothesis = h, .material = m});
}  // end of Behaviour_loadFromDatabase

void declareModel(pybind11::module_& m) {
  using namespace pybind11::literals;
  // wrapping free functions
  m.def("load", pybind11::overload_cast<const std::string&, const std::string&,
                                        const mgis::behaviour::Hypothesis>(
                    mgis::model::load));
  m.def("loadFromDatabase", Model_loadFromDatabase, "name"_a, "hypothesis"_a,
        "material"_a = std::optional<std::string>{});

}  // end of declareModel
