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
     * \param[in] i: cell index
     */
    real operator()(const size_type) const;

   private:
    //! \brief evaluator of deformation gradient
    DeformationGradientEvaluatorType deformation_gradient_evaluator;
    //! \brief evaluator of the first Piola-Kirchhoff stress
    PK1EvaluatorType pk1_evaluator;
  };

}  // end of namespace mgis::function

#include "MGIS/Function/MechanicalEvaluators/CauchyStressFromFirstPiolaKirchhoffStressEvaluator.ixx"

#endif /* LIB_CAUCHYSTRESSFROMFIRSTPIOLAKIRCHHOFFSTRESSEVALUATOR_HXX */
