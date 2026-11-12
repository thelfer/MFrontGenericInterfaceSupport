/*!
 * \file   MGIS/Function/EvaluatorModifierBase.hxx
 * \brief File which defines the `EvaluatorModifierBase` class, a helper base
 * class used to build a new evaluator by applying a modifier to the values
 * returned by an existing evaluator.
 * \author Thomas Helfer
 * \date   07/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_EVALUATORMODIFIERBASE_HXX
#define LIB_MGIS_FUNCTION_EVALUATORMODIFIERBASE_HXX

#include "MGIS/AbstractErrorHandler.hxx"
#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/EvaluatorConcept.hxx"
#include "MGIS/Function/EvaluatorModifierConcept.hxx"

namespace mgis::function {

  /*!
   * \brief a base class to construct a new evaluator by applying a modifier to
   * the values returned by an existing evaluator
   *
   * This class implements the CRTP idiom: the call operators forward the
   * values computed by the underlying evaluator to the `apply` member
   * function of the derived class, which is responsible for modifying them.
   * Only the call operators which match the way the underlying evaluator can
   * be called are available.
   *
   * \tparam Child: type of the derived class, which must implement an `apply`
   * member function accepting the values returned by the underlying evaluator
   * \tparam EvaluatorType: type of the modified evaluator
   */
  template <typename Child, EvaluatorConcept EvaluatorType>
  struct EvaluatorModifierBase {
    //! \brief type of the space on which the underlying evaluator is defined
    using Space = evaluator_space<EvaluatorType>;
    /*!
     * \brief constructor
     * \param[in] e: modified evaluator
     */
    constexpr EvaluatorModifierBase(const EvaluatorType&);
    //! \brief copy constructor
    constexpr EvaluatorModifierBase(const EvaluatorModifierBase&);
    //! \brief move constructor
    constexpr EvaluatorModifierBase(EvaluatorModifierBase&&);
    /*!
     * \brief return the underlying space
     * \return the space on which the underlying evaluator is defined
     */
    [[nodiscard]] constexpr decltype(auto) getSpace() const;
    /*!
     * \brief perform consistency checks
     * \param[in] ctx: error handler used to report the errors
     * \return true if the values computed by the underlying evaluator are valid
     */
    [[nodiscard]] constexpr bool check(AbstractErrorHandler&) const;
    /*!
     * \brief call operator: the value computed by the underlying evaluator is
     * passed to the `apply` member function of the derived class
     * \param[in] e: element index, i.e. index of the integration point
     * \return the modified value
     * \note only available if the underlying evaluator can be called with an
     * element index
     */
    [[nodiscard]] constexpr auto operator()(const element_index<Space>&) const
        requires(internals::EvaluatorResultQuery<EvaluatorType>::b1);
    /*!
     * \brief call operator: the value computed by the underlying evaluator is
     * passed to the `apply` member function of the derived class
     * \param[in] wk: element workspace
     * \param[in] e: element index, i.e. index of the integration point
     * \return the modified value
     * \note only available if the underlying evaluator can be called with an
     * element workspace and an element index
     */
    [[nodiscard]] constexpr auto operator()(const element_workspace<Space>&,
                                            const element_index<Space>&) const
        requires(internals::EvaluatorResultQuery<EvaluatorType>::b2);
    /*!
     * \brief call operator: the value computed by the underlying evaluator is
     * passed to the `apply` member function of the derived class
     * \param[in] e: cell index
     * \param[in] i: quadrature point index
     * \return the modified value
     * \note only available if the underlying evaluator can be called with a
     * cell index and a quadrature point index
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_index<Space>&, const quadrature_point_index<Space>&) const
        requires(internals::EvaluatorResultQuery<EvaluatorType>::b3);
    /*!
     * \brief call operator: the value computed by the underlying evaluator is
     * passed to the `apply` member function of the derived class
     * \param[in] wk: cell workspace
     * \param[in] e: cell index
     * \param[in] i: quadrature point index
     * \return the modified value
     * \note only available if the underlying evaluator can be called with a
     * cell workspace, a cell index and a quadrature point index
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_workspace<Space>&,
        const cell_index<Space>&,
        const quadrature_point_index<Space>&) const
        requires(internals::EvaluatorResultQuery<EvaluatorType>::b4);

   protected:
    //! \brief underlying evaluator
    EvaluatorType evaluator;
  };

  /*!
   * \brief return the underlying space
   * \tparam Child: type of the derived class
   * \tparam EvaluatorType: type of the modified evaluator
   * \param[in] e: evaluator modifier
   * \return the space on which the modifier is defined
   */
  template <typename Child, EvaluatorConcept EvaluatorType>
  [[nodiscard]] constexpr decltype(auto) getSpace(
      const EvaluatorModifierBase<Child, EvaluatorType>&);
  /*!
   * \brief perform consistency checks
   * \tparam Child: type of the derived class
   * \tparam EvaluatorType: type of the modified evaluator
   * \param[in] eh: error handler used to report the errors
   * \param[in] e: evaluator modifier
   * \return true if the values computed by the modifier are valid
   */
  template <typename Child, EvaluatorConcept EvaluatorType>
  [[nodiscard]] constexpr bool check(
      AbstractErrorHandler&,
      const EvaluatorModifierBase<Child, EvaluatorType>&);

}  // namespace mgis::function

#include "MGIS/Function/EvaluatorModifierBase.ixx"

#endif /* LIB_MGIS_FUNCTION_EVALUATORMODIFIERBASE_HXX */
