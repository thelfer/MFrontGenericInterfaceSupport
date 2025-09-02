/*!
 * \file   MGIS/Function/EvaluatorModifierBase.hxx
 * \brief
 * \author Thomas Helfer
 * \date   07/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_EVALUATORMODIFIERBASE_HXX
#define LIB_MGIS_FUNCTION_EVALUATORMODIFIERBASE_HXX

#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/Evaluator.hxx"

namespace mgis::function {

  /*!
   * \brief a base class to construct a new evaluator by applying a modifier to
   * the values returned by an existing evaluator
   * \tparam Child: child class
   * \tparam EvaluatorType: modified evaluator
   */
  template <typename Child, EvaluatorConcept EvaluatorType>
  struct EvaluatorModifierBase {
    //! \brief a simple alias
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
        requires(internals::EvaluatorResultQuery<EvaluatorType>::b1);
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    constexpr auto operator()(const element_workspace<Space>&,
                              const element_index<Space>&) const
        requires(internals::EvaluatorResultQuery<EvaluatorType>::b2);
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    constexpr auto operator()(const cell_index<Space>&,
                              const quadrature_point_index<Space>&) const
        requires(internals::EvaluatorResultQuery<EvaluatorType>::b3);
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    constexpr auto operator()(const cell_workspace<Space>&,
                              const cell_index<Space>&,
                              const quadrature_point_index<Space>&) const
        requires(internals::EvaluatorResultQuery<EvaluatorType>::b4);

   protected:
    //! \brief underlying evaluator
    EvaluatorType evaluator;
  };

  template <typename Child, EvaluatorConcept EvaluatorType>
  constexpr decltype(auto) getSpace(
      const EvaluatorModifierBase<Child, EvaluatorType>&);
  //! \brief allocate internal workspace
  template <typename Child, EvaluatorConcept EvaluatorType>
  constexpr void allocateWorkspace(
      EvaluatorModifierBase<Child, EvaluatorType>&);

}  // namespace mgis::function

#include "MGIS/Function/EvaluatorModifierBase.ixx"

#endif /* LIB_MGIS_FUNCTION_EVALUATORMODIFIERBASE_HXX */
