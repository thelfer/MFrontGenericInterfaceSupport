/*!
 * \file   include/MGIS/Model/Model.hxx
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

#ifndef LIB_MGIS_MODEL_MODEL_HXX
#define LIB_MGIS_MODEL_MODEL_HXX

#include "MGIS/Config.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"

namespace mgis::model {

  /*!
   * \brief a simple alias
   *
   * In `MGIS`, a model is a behaviour without gradients, thermodynamic forces
   * nor tangent operator blocks.
   */
  using Model = mgis::behaviour::Behaviour;

  /*!
   * \brief load the description of a model from a library
   *
   * \param[in] l: library name
   * \param[in] m: model name
   * \param[in] h: modelling hypothesis
   * \return the model description
   * \note: use of `std::string` rather than `mgis::string_view` is
   * meaningfull here
   */
  MGIS_EXPORT Model load(const std::string &,
                         const std::string &,
                         const mgis::behaviour::Hypothesis);
  /*!
   * \brief load the description of a model from a library
   *
   * \param[in] l: library name
   * \param[in] m: model name
   * \param[in] h: modelling hypothesis
   * \return the model description
   * \note: use of `std::string` rather than `mgis::string_view` is
   * meaningfull here
   */
  MGIS_EXPORT std::optional<Model> load(
      Context &,
      const std::string &,
      const std::string &,
      const mgis::behaviour::Hypothesis) noexcept;
  /*!
   * \brief load a model from the database
   *
   * \return the model description
   */
  MGIS_EXPORT Model
  loadFromDatabase(const mgis::behaviour::LoadFromDatabaseOptions &);
  /*!
   * \brief load a model from the database
   *
   * \param[in] ctx: execution context
   * \param[in] opts: options to select the model
   * \return the model description
   */
  MGIS_EXPORT std::optional<Model> loadFromDatabase(
      Context &, const mgis::behaviour::LoadFromDatabaseOptions &) noexcept;

}  // namespace mgis::model

#endif /* LIB_MGIS_MODEL_MODEL_HXX */
