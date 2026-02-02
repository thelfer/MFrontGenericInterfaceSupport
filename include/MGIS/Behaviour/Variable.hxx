/*!
 * \file   MGIS/Behaviour/Variable.hxx
 * \brief
 * \author Thomas Helfer
 * \date   19/06/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_BEHAVIOUR_VARIABLE_HXX
#define LIB_MGIS_BEHAVIOUR_VARIABLE_HXX

#include <string>
#include <vector>
#include <string_view>
#include "MGIS/Config.hxx"
#include "MGIS/Context.hxx"
#include "MGIS/Behaviour/Hypothesis.hxx"

namespace mgis::behaviour {

  /*!
   * \brief structure describing a variable
   */
  struct Variable {
    //! \brief name of the variable
    std::string name;
    //! \brief type of the variable
    enum Type {
      SCALAR = 0,
      VECTOR = 2,
      VECTOR_1D = 10,
      VECTOR_2D = 18,
      VECTOR_3D = 26,
      STENSOR = 1,
      STENSOR_1D = 9,
      STENSOR_2D = 17,
      STENSOR_3D = 25,
      TENSOR = 3,
      TENSOR_1D = 11,
      TENSOR_2D = 19,
      TENSOR_3D = 27,
      HIGHER_ORDER_TENSOR = 4,
      ARRAY = 5
    } type;
    //! brief type identifier
    int type_identifier = 0;
  };  // end of struct Variable

  /*!
   * \return a boolean stating that a variable with the given name
   * is in the container.
   * \param[in] vs: variables
   * \param[in] n: name
   */
  MGIS_EXPORT [[nodiscard]] bool contains(const std::vector<Variable> &,
                                          const std::string_view) noexcept;
  /*!
   * \return the variable with the given name
   * \param[in] vs: variables
   * \param[in] n: name
   */
  MGIS_EXPORT [[nodiscard]] const Variable &getVariable(
      const std::vector<Variable> &, const std::string_view);
  /*!
   * \return the variable with the given name
   * \param[in] vs: variables
   * \param[in] n: name
   */
  MGIS_EXPORT [[nodiscard]] std::optional<const Variable *> getVariable(
      Context &,
      const std::vector<Variable> &,
      const std::string_view) noexcept;
  /*!
   * \return the type of a variable from an identifier
   * \param[in] id: type identifier
   */
  MGIS_EXPORT [[nodiscard]] Variable::Type getVariableType(const int);
  /*!
   * \return a symbolic representation from a type identifier
   * \param[in] id: type identifier
   */
  MGIS_EXPORT [[nodiscard]] std::string getVariableTypeSymbolicRepresentation(
      const int);
  /*!
   * \return the size of a variable
   * \param[in] v: variable
   * \param[in] h: modelling hypothesis
   */
  MGIS_EXPORT [[nodiscard]] size_type getVariableSize(const Variable &,
                                                      const Hypothesis);
  /*!
   * \return the size of a variable
   * \param[in, out] ctx: execution context
   * \param[in] v: variable
   * \param[in] h: modelling hypothesis
   */
  MGIS_EXPORT [[nodiscard]] std::optional<size_type> getVariableSize(
      Context &, const Variable &, const Hypothesis) noexcept;
  /*!
   * \return the size of an array that may contain the values described by the
   * given array of variables
   * \param[in] vs: variables
   * \param[in] h: modelling hypothesis
   */
  MGIS_EXPORT [[nodiscard]] size_type getArraySize(
      const std::vector<Variable> &, const Hypothesis);
  /*!
   * \return the offset of the given variable for the given hypothesis
   * \param[in] vs: variables
   * \param[in] n: variable name
   * \param[in] h: modelling hypothesis
   */
  MGIS_EXPORT [[nodiscard]] size_type getVariableOffset(
      const std::vector<Variable> &, const std::string_view, const Hypothesis);
  /*!
   * \return the type of the given variable as a string
   * \param[in] v: variable
   */
  MGIS_EXPORT [[nodiscard]] std::string getVariableTypeAsString(
      const Variable &);

}  // end of namespace mgis::behaviour

#endif /* LIB_MGIS_BEHAVIOUR_VARIABLE_HXX */
