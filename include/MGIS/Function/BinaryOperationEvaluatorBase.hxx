/*!
 * \file   MGIS/Function/BinaryOperationEvaluatorBase.hxx
 * \brief
 * \author Thomas Helfer
 * \date   07/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATORBASE_HXX
#define LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATORBASE_HXX

#include "MGIS/AbstractErrorHandler.hxx"
#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/Evaluator.hxx"

namespace mgis::function {

  /*!
   * \brief the set of requirements that must be satisfied by both operands of
   * a binary operation on evaluators
   */
  template <EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  inline constexpr auto BinaryOperationEvaluatorBaseRequirement =
      (internals::same_decay_type<
          decltype(getSpace(std::declval<FirstEvaluatorType>())),
          decltype(getSpace(std::declval<FirstEvaluatorType>()))>)&&  //
      (((ElementEvaluatorConcept<FirstEvaluatorType>)&&(
           ElementEvaluatorConcept<SecondEvaluatorType>)) ||
       ((QuadratureEvaluatorConcept<FirstEvaluatorType>)&&(
           QuadratureEvaluatorConcept<SecondEvaluatorType>)));

  /*!
   * \brief a base class to construct an evaluator by applying a binary
   * operation on the results of two evaluators.
   * \tparam Child: child class
   * \tparam FirstEvaluatorType: evaluator used to compute
   * the first argument of the binary operation
   * \tparam SecondEvaluatorType: evaluator used to compute
   * the first second of the binary operation
   */
  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  requires(BinaryOperationEvaluatorBaseRequirement<FirstEvaluatorType,
                                                   SecondEvaluatorType>)  //
      struct BinaryOperationEvaluatorBase {
    //! \brief a simple alias
    using Space = evaluator_space<FirstEvaluatorType>;
    // boolean stating if both evaluators matches the ElementEvaluatorConcept
    static constexpr auto isElementEvaluator =
        (ElementEvaluatorConcept<FirstEvaluatorType>)&&(
            ElementEvaluatorConcept<SecondEvaluatorType>);
    // boolean stating if both evaluators matches the QuadratureEvaluatorConcept
    static constexpr auto isQuadratureEvaluator =
        (QuadratureEvaluatorConcept<FirstEvaluatorType>)&&(
            QuadratureEvaluatorConcept<SecondEvaluatorType>);
    /*!
     * \brief constructor
     * \param[in] e1: first evaluator
     * \param[in] e2: second evaluator
     */
    constexpr BinaryOperationEvaluatorBase(const FirstEvaluatorType&,
                                           const SecondEvaluatorType&);
    //! \brief copy constructor
    constexpr BinaryOperationEvaluatorBase(const BinaryOperationEvaluatorBase&);
    //! \brief move constructor
    constexpr BinaryOperationEvaluatorBase(BinaryOperationEvaluatorBase&&);
    //! \brief perform consistency checks
    constexpr bool check(AbstractErrorHandler&) const;
    //! \brief allocate internal workspace
    constexpr void allocateWorkspace();
    //! \brief return the underlying space
    constexpr decltype(auto) getSpace() const;
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    constexpr auto operator()(const element_index<Space>&) const
        requires((internals::EvaluatorResultQuery<FirstEvaluatorType>::b1) &&
                 (internals::EvaluatorResultQuery<SecondEvaluatorType>::b1));
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    constexpr auto operator()(const element_workspace<Space>&,
                              const element_index<Space>&) const
        requires((internals::EvaluatorResultQuery<FirstEvaluatorType>::b2) &&
                 (internals::EvaluatorResultQuery<SecondEvaluatorType>::b2));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    constexpr auto operator()(const cell_index<Space>,
                              const quadrature_point_index<Space>) const
        requires((internals::EvaluatorResultQuery<FirstEvaluatorType>::b3) &&
                 (internals::EvaluatorResultQuery<SecondEvaluatorType>::b3));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    constexpr auto operator()(const cell_workspace<Space>&,
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

  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  constexpr decltype(auto) getSpace(
      const BinaryOperationEvaluatorBase<Child,
                                         FirstEvaluatorType,
                                         SecondEvaluatorType>&);
  //! \brief allocate internal workspace
  template <typename Child,
            EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  constexpr void allocateWorkspace(
      BinaryOperationEvaluatorBase<Child,
                                   FirstEvaluatorType,
                                   SecondEvaluatorType>&);

}  // end of namespace mgis::function

#include "MGIS/Function/BinaryOperationEvaluatorBase.ixx"

#endif /* LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATORBASE_HXX */
