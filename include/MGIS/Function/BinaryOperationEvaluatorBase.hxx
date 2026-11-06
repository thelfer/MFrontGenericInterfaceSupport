/*!
 * \file   MGIS/Function/BinaryOperationEvaluatorBase.hxx
 * \brief  This file declares the `BinaryOperationEvaluatorBase` class,
 * a base class used to build an evaluator by applying a binary operation
 * to the values returned by two evaluators.
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

#ifndef LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATORBASE_HXX
#define LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATORBASE_HXX

#include "MGIS/AbstractErrorHandler.hxx"
#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/EvaluatorConcept.hxx"

namespace mgis::function {

  /*!
   * \brief the set of requirements that must be satisfied by both operands of
   * a binary operation on evaluators
   *
   * The two evaluators must be defined on the same space and must either both
   * match the `ElementEvaluatorConcept` concept or both the
   * `QuadratureEvaluatorConcept` concept: the values of the two operands are
   * always computed at the same integration point.
   *
   * \tparam FirstEvaluatorType: evaluator used to compute the first argument
   * of the binary operation
   * \tparam SecondEvaluatorType: evaluator used to compute the second
   * argument of the binary operation
   */
  template <EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  inline constexpr auto BinaryOperationEvaluatorBaseRequirement =
      (std::same_as<evaluator_space<FirstEvaluatorType>,
                    evaluator_space<SecondEvaluatorType>>)&&  //
      (((ElementEvaluatorConcept<FirstEvaluatorType>)&&(
           ElementEvaluatorConcept<SecondEvaluatorType>)) ||
       ((QuadratureEvaluatorConcept<FirstEvaluatorType>)&&(
           QuadratureEvaluatorConcept<SecondEvaluatorType>)));

  /*!
   * \brief a base class to construct an evaluator by applying a binary
   * operation on the results of two evaluators
   *
   * This class is meant to be used with the curiously recurring template
   * pattern: the child class gives the binary operation through an `apply`
   * member function taking the values returned by the two evaluators and
   * returning the values of the resulting evaluator.
   *
   * The call operators of this class forward the values computed by the two
   * evaluators to the `apply` member function of the child class: only the
   * call operators which match the way both evaluators can be called are
   * available.
   *
   * \tparam Child: type of the class deriving from this base class
   * \tparam FirstEvaluatorType: evaluator used to compute the first argument
   * of the binary operation
   * \tparam SecondEvaluatorType: evaluator used to compute the second
   * argument of the binary operation
   * \note both evaluators must satisfy the requirements stated by
   * `BinaryOperationEvaluatorBaseRequirement`
   */
  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationEvaluatorBaseRequirement<FirstEvaluatorType,
                                                   SecondEvaluatorType>)  //
      struct BinaryOperationEvaluatorBase {
    //! \brief type of the space on which the two evaluators are defined
    using Space = evaluator_space<FirstEvaluatorType>;
    //! \brief true if both evaluators match the `ElementEvaluatorConcept`
    static constexpr auto isElementEvaluator =
        (ElementEvaluatorConcept<FirstEvaluatorType>)&&(
            ElementEvaluatorConcept<SecondEvaluatorType>);
    //! \brief true if both evaluators match the `QuadratureEvaluatorConcept`
    static constexpr auto isQuadratureEvaluator =
        (QuadratureEvaluatorConcept<FirstEvaluatorType>)&&(
            QuadratureEvaluatorConcept<SecondEvaluatorType>);
    /*!
     * \brief constructor
     * \param[in] e1: evaluator used to compute the first argument of the
     * binary operation, a copy of it is stored in the object
     * \param[in] e2: evaluator used to compute the second argument of the
     * binary operation, a copy of it is stored in the object
     */
    constexpr BinaryOperationEvaluatorBase(const FirstEvaluatorType&,
                                           const SecondEvaluatorType&);
    //! \brief copy constructor
    constexpr BinaryOperationEvaluatorBase(const BinaryOperationEvaluatorBase&);
    //! \brief move constructor
    constexpr BinaryOperationEvaluatorBase(BinaryOperationEvaluatorBase&&);
    /*!
     * \brief return the underlying space
     * \return the space on which the two evaluators are defined
     */
    [[nodiscard]] constexpr decltype(auto) getSpace() const;
    /*!
     * \brief perform consistency checks
     * \param[in] ctx: error handler used to report the errors
     * \return true if the values computed by the two evaluators are valid,
     * the errors being reported using the error handler `ctx`
     */
    [[nodiscard]] constexpr bool check(AbstractErrorHandler&) const;
    /*!
     * \brief call operator
     * \param[in] i: index of the element on which the values are computed
     * \return values computed by the `apply` member function of the child class
     * from the values returned by the two evaluators at the integration point
     * `i`, which is the index of the element on which the values are computed
     * \note only available if both evaluators can be called with an element
     * index
     */
    [[nodiscard]] constexpr auto operator()(const element_index<Space>&) const
        requires((internals::EvaluatorResultQuery<FirstEvaluatorType>::b1) &&
                 (internals::EvaluatorResultQuery<SecondEvaluatorType>::b1));
    /*!
     * \brief call operator
     * \param[in] wk: workspace of the element on which the values are computed
     * \param[in] i: index of the integration point of this element
     * \return values computed by the `apply` member function of the child class
     * from the values returned by the two evaluators at the integration point
     * `i` of the element whose workspace is `wk`
     * \note only available if both evaluators can be called with an element
     * workspace and an element index
     */
    [[nodiscard]] constexpr auto operator()(const element_workspace<Space>&,
                                            const element_index<Space>&) const
        requires((internals::EvaluatorResultQuery<FirstEvaluatorType>::b2) &&
                 (internals::EvaluatorResultQuery<SecondEvaluatorType>::b2));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: quadrature point index of the cell of index `e`
     * \return values computed by the `apply` member function of the child class
     * from the values returned by the two evaluators at the quadrature point
     * `i` of the cell `e`
     * \note only available if both evaluators can be called with a cell index
     * and a quadrature point index
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_index<Space>, const quadrature_point_index<Space>) const
        requires((internals::EvaluatorResultQuery<FirstEvaluatorType>::b3) &&
                 (internals::EvaluatorResultQuery<SecondEvaluatorType>::b3));
    /*!
     * \brief call operator
     * \param[in] wk: workspace of the cell of index `e`
     * \param[in] e: cell index
     * \param[in] i: quadrature point index of the cell of index `e`
     * \return values computed by the `apply` member function of the child class
     * from the values returned by the two evaluators at the quadrature point
     * `i` of the cell `e`, whose workspace is `wk`
     * \note only available if both evaluators can be called with a cell
     * workspace, a cell index and a quadrature point index
     */
    [[nodiscard]] constexpr auto operator()(
        const cell_workspace<Space>&,
        const cell_index<Space>,
        const quadrature_point_index<Space>) const
        requires((internals::EvaluatorResultQuery<FirstEvaluatorType>::b4) &&
                 (internals::EvaluatorResultQuery<SecondEvaluatorType>::b4));

   protected:
    //! \brief evaluator of the first argument of the binary operation
    FirstEvaluatorType first_evaluator;
    //! \brief evaluator of the second argument of the binary operation
    SecondEvaluatorType second_evaluator;
  };

  /*!
   * \brief return the underlying space
   * \tparam Child: type of the class deriving from this base class
   * \tparam FirstEvaluatorType: evaluator used to compute the first
   * argument of the binary operation
   * \tparam SecondEvaluatorType: evaluator used to compute the second
   * argument of the binary operation
   * \param[in] e: evaluator of the result of the binary operation
   * \return the space on which the two evaluators are defined
   */
  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  [[nodiscard]] constexpr decltype(auto) getSpace(
      const BinaryOperationEvaluatorBase<Child,
                                         FirstEvaluatorType,
                                         SecondEvaluatorType>&);
  /*!
   * \brief perform consistency checks
   * \tparam Child: type of the class deriving from this base class
   * \tparam FirstEvaluatorType: evaluator used to compute the first
   * argument of the binary operation
   * \tparam SecondEvaluatorType: evaluator used to compute the second
   * argument of the binary operation
   * \param[in] eh: error handler used to report the errors
   * \param[in] e: evaluator of the result of the binary operation
   * \return true if the values computed by the evaluator are valid
   */
  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  [[nodiscard]] constexpr bool check(
      AbstractErrorHandler&,
      const BinaryOperationEvaluatorBase<Child,
                                         FirstEvaluatorType,
                                         SecondEvaluatorType>&);

}  // end of namespace mgis::function

#include "MGIS/Function/BinaryOperationEvaluatorBase.ixx"

#endif /* LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATORBASE_HXX */
