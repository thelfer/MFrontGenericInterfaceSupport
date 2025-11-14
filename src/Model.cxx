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
#include "MGIS/Context.hxx"
#include "MGIS/Model/Model.hxx"

namespace mgis::model {

  static void checkModel(const Model &m) {
    auto throw_if = [&m](const bool c, const char *const type) {
      if (c) {
        auto msg = std::string{
            "mgis::model::loadModel: "
            "model '" +
            m.function + "' in library '" + m.library +
            "' shall not declare any "};
        msg += type;
        mgis::raise(type);
      }
    };
    throw_if(!m.gradients.empty(), "gradient");
    throw_if(!m.thermodynamic_forces.empty(), "thermodynamic force");
    throw_if(!m.to_blocks.empty(), "tangent operator block");
  }

  Model load(const std::string &l,
             const std::string &m,
             const mgis::behaviour::Hypothesis h) {
    const auto model = ::mgis::behaviour::load(l, m, h);
    checkModel(model);
    return model;
  }  // end of load

  std::optional<Model> load(Context &ctx,
                            const std::string &l,
                            const std::string &m,
                            const mgis::behaviour::Hypothesis h) noexcept {
    try {
      return ::mgis::model::load(l, m, h);
    } catch (...) {
      registerExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of load

  Model loadFromDatabase(const mgis::behaviour::LoadFromDatabaseOptions &opts) {
    const auto model = mgis::behaviour::loadFromDatabase(opts);
    checkModel(model);
    return model;
  }  // end of load

  std::optional<Model> load(
      Context &ctx,
      const mgis::behaviour::LoadFromDatabaseOptions &opts) noexcept {
    try {
      return ::mgis::model::loadFromDatabase(opts);
    } catch (...) {
      registerExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of load

}  // namespace mgis::model
