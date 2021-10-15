/*!
 * \file   src/Model.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/10/2021
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include "MGIS/Raise.hxx"
#include "MGIS/Model/Model.hxx"

namespace mgis::model {

  Model load(const std::string &l,
             const std::string &m,
             const mgis::behaviour::Hypothesis h) {
    const auto model = mgis::behaviour::load(l, m, h);
    auto throw_if = [&l, &m](const bool c, const char *const type) {
      if (c) {
        auto msg = std::string{
            "mgis::model::loadModel: "
            "model '" +
            m + "' in library '" + l + "' shall not declare any "};
        msg += type;
        mgis::raise(type);
      }
    };
    throw_if(!model.gradients.empty(), "gradient");
    throw_if(!model.thermodynamic_forces.empty(), "thermodynamic force");
    throw_if(!model.to_blocks.empty(), "tangent operator block");
    return model;
  }  // end of load

}  // namespace mgis::model
