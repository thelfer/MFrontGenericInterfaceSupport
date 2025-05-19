/*!
 * \file   MGIS/Function/BinaryOperationEvaluatorBase.hxx
 * \brief
 * \author Thomas Helfer
 * \date   07/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATORBASE_HXX
#define LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATORBASE_HXX

#include "MGIS/Context.hxx"
#include "MGIS/Function/Space.hxx"
#include "MGIS/Function/Evaluator.hxx"

namespace mgis::function {

  template <EvaluatorConcept FirstEvaluatorType,
            EvaluatorConcept SecondEvaluatorType>
  inline constexpr auto BinaryOperationEvaluatorBaseRequirement =
      (internals::same_decay_type<
          decltype(std::declval<FirstEvaluatorType>().getSpace()),
          decltype(std::declval<FirstEvaluatorType>().getSpace())>)&&  //
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
    using Space =
        std::decay_t<decltype(std::declval<FirstEvaluatorType>().getSpace())>;
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
    BinaryOperationEvaluatorBase(const FirstEvaluatorType&,
                                 const SecondEvaluatorType&);
    //! \brief copy constructor
    BinaryOperationEvaluatorBase(const BinaryOperationEvaluatorBase&);
    //! \brief move constructor
    BinaryOperationEvaluatorBase(BinaryOperationEvaluatorBase&&);
    //! \brief perform consistency checks
    bool check(Context&) const;
    //! \brief allocate internal workspace
    void allocateWorkspace();
    //! \brief return the underlying space
    const auto& getSpace() const;
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    auto operator()(const element_index<Space>&) const
        requires((internals::EvaluatorResultQuery<FirstEvaluatorType>::b1) &&
                 (internals::EvaluatorResultQuery<SecondEvaluatorType>::b1));
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    auto operator()(const element_workspace<Space>&,
                    const element_index<Space>&) const
        requires((internals::EvaluatorResultQuery<FirstEvaluatorType>::b2) &&
                 (internals::EvaluatorResultQuery<SecondEvaluatorType>::b2));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    auto operator()(const cell_index<Space>,
                    const quadrature_point_index<Space>) const
        requires((internals::EvaluatorResultQuery<FirstEvaluatorType>::b3) &&
                 (internals::EvaluatorResultQuery<SecondEvaluatorType>::b3));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    auto operator()(const cell_workspace<Space>&,
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

}  // end of namespace mgis::function

#include "MGIS/Function/BinaryOperationEvaluatorBase.ixx"

#endif /* LIB_MGIS_FUNCTION_BINARYOPERATIONEVALUATORBASE_HXX */
