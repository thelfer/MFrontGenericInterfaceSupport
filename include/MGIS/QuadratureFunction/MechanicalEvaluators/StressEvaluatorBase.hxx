/*!
 * \file   StressEvaluatorBase.hxx
 * \brief
 * \author Thomas Helfer
 * \date   02/05/2025
 */

#ifndef MGIS_HAVE_TFEL
#error "TFEL is required to use mechanical evaluators"
#endif /* MGIS_HAVE_TFEL */

#ifndef LIB_MGIS_QUADRATUREFUNCTION_STRESSEVALUATORBASE_HXX
#define LIB_MGIS_QUADRATUREFUNCTION_STRESSEVALUATORBASE_HXX

#include "MGIS/QuadratureFunction/Evaluators.hxx"

namespace mgis::quadrature_function {

  /*!
   * \brief a base class for evaluators modifying a stress tensor
   * \tparam N: space dimension
   * \tparam StressEvaluatorType: evaluator of the stress
   * \tparam symmetric: boolean stating of the stress tensor is symmetric
   */
  template <unsigned short N,
            EvaluatorConcept StressEvaluatorType,
            bool isSymmetric>
  struct StressEvaluatorBase {
    /*!
     * \brief constructor
     * \param[in] e: evaluator of the stress
     */
    StressEvaluatorBase(const StressEvaluatorType&);
    //! \brief copy constructor
    StressEvaluatorBase(const StressEvaluatorBase&);
    //! \brief move constructor
    StressEvaluatorBase(StressEvaluatorBase&&);
    //! \brief perform consistency checks
    void check() const;
    //! \brief allocate internal workspace
    void allocateWorkspace();
    //! \brief return the underlying partial quadrature space
    const AbstractQuadratureSpace& getQuadratureSpace() const;
    //! \return the number of components
    constexpr size_type getNumberOfComponents() const noexcept;

   private:
    //! \brief evaluator of the stress
    StressEvaluatorType stress_evaluator;
  };

}  // namespace mgis::quadrature_function

#include "MGIS/QuadratureFunction/MechanicalEvaluators/StressEvaluatorBase.ixx"

#endif /* LIB_MGIS_QUADRATUREFUNCTION_STRESSEVALUATORBASE_HXX */
