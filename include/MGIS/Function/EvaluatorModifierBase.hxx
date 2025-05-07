/*!
 * \file   MGIS/Function/EvaluatorModifierBase.hxx
 * \brief
 * \author Thomas Helfer
 * \date   07/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_EVALUATORMODIFIERBASE_HXX
#define LIB_MGIS_FUNCTION_EVALUATORMODIFIERBASE_HXX

#include "MGIS/Function/Space.hxx"
#include "MGIS/Function/Evaluators.hxx"

namespace mgis::function {

  /*!
   * \brief a base class for applying a modifier to an evaluator
   * \tparam Child: child class
   * \tparam EvaluatorType: evaluator of the stress
   */
  template <typename Child, EvaluatorConcept EvaluatorType>
  struct EvaluatorModifierBase {
    //! \brief a simple alias
    using Space =
        std::decay_t<decltype(std::declval<EvaluatorType>().getSpace())>;
    /*!
     * \brief constructor
     * \param[in] e: evaluator of the stress
     */
    EvaluatorModifierBase(const EvaluatorType&);
    //! \brief copy constructor
    EvaluatorModifierBase(const EvaluatorModifierBase&);
    //! \brief move constructor
    EvaluatorModifierBase(EvaluatorModifierBase&&);
    //! \brief perform consistency checks
    void check() const;
    //! \brief allocate internal workspace
    void allocateWorkspace();
    //! \brief return the underlying partial quadrature space
    const auto& getSpace() const;
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    auto operator()(const element_index<Space>&) const
        requires(ElementSpaceConcept<Space> && !(hasElementWorkspace<Space>));
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    auto operator()(const element_workspace<Space>&,
                    const element_index<Space>&) const
        requires(ElementSpaceConcept<Space>&& hasElementWorkspace<Space>);
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    auto operator()(const cell_index<Space>,
                    const quadrature_point_index<Space>) const
        requires(QuadratureSpaceConcept<Space> && (!hasCellWorkspace<Space>));
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    auto operator()(const cell_workspace<Space>&,
                    const cell_index<Space>,
                    const quadrature_point_index<Space>) const
        requires(QuadratureSpaceConcept<Space>&& hasCellWorkspace<Space>);

   protected:
    //! \brief underlying evaluator
    EvaluatorType evaluator;
  };

}  // namespace mgis::function

#include "MGIS/Function/EvaluatorModifierBase.ixx"

#endif /* LIB_MGIS_FUNCTION_EVALUATORMODIFIERBASE_HXX */
