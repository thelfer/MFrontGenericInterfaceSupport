/*!
 * \file   CauchyStressFromFirstPiolaKirchhoffStressEvaluator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#ifndef MGIS_HAVE_TFEL
#error "TFEL is required to use mechanical evaluators"
#endif /* MGIS_HAVE_TFEL */

#ifndef LIB_MGIS_CAUCHYSTRESSFROMFIRSTPIOLAKIRCHHOFFSTRESSEVALUATOR_HXX
#define LIB_MGIS_CAUCHYSTRESSFROMFIRSTPIOLAKIRCHHOFFSTRESSEVALUATOR_HXX

#include "MGIS/Function/Evaluators.hxx"

namespace mgis::function {

  /*!
   * \brief an evaluator returning the Cauchy stress stress using evaluators of
   * the deformation gradient and the first Piola-Kirchhoff stress
   *
   * \tparam N: space dimension
   * \tparam DeformationGradientEvaluatorType: evaluator of the deformation
   * gradient
   * \tparam PK1EvaluatorType: evaluator of the first Piola-Kirchhoff stress
   */
  template <unsigned short N,
            EvaluatorConcept DeformationGradientEvaluatorType,
            EvaluatorConcept PK1EvaluatorType>
  struct CauchyStressFromFirstPiolaKirchhoffStressEvaluator {
    //! \brief a simple alias
    using Space = std::decay_t<
        decltype(std::declval<std::decay_t<PK1EvaluatorType>>().getSpace())>;
    /*!
     * \brief constructor
     * \param[in] e1: evaluator of the deformation gradient
     * \param[in] e2: evaluator of the first Piola-Kirchhoff stress
     */
    CauchyStressFromFirstPiolaKirchhoffStressEvaluator(
        const DeformationGradientEvaluatorType&, const PK1EvaluatorType&);
    //! \brief perform consistency checks
    void check() const;
    //! \brief allocate internal workspace
    void allocateWorkspace();
    //! \brief return the underlying partial quadrature space
    const auto& getSpace() const;
    //! \return the number of components
    constexpr size_type getNumberOfComponents() const noexcept;
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

   private:
    //
    auto apply(const auto& pk1_values, const auto& F_values) const;
    //! \brief evaluator of deformation gradient
    DeformationGradientEvaluatorType deformation_gradient_evaluator;
    //! \brief evaluator of the first Piola-Kirchhoff stress
    PK1EvaluatorType pk1_evaluator;
  };

}  // end of namespace mgis::function

#include "MGIS/Function/MechanicalEvaluators/CauchyStressFromFirstPiolaKirchhoffStressEvaluator.ixx"

#endif /* LIB_CAUCHYSTRESSFROMFIRSTPIOLAKIRCHHOFFSTRESSEVALUATOR_HXX */
