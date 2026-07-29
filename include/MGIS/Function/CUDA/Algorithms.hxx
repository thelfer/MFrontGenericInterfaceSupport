/*!
 * \file   MGIS/Function/CUDA/Algorithms.hxx
 * \brief
 * \author Thomas Helfer
 * \date   28/06/2026
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_CUDA_ALGORITHMS_HXX
#define LIB_MGIS_FUNCTION_CUDA_ALGORITHMS_HXX

#include "MGIS/Config.hxx"
#include "MGIS/Context.hxx"
#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/EvaluatorConcept.hxx"
#include "MGIS/Function/FunctionConcept.hxx"

namespace mgis::function {

  /*!
   * \brief execution configuration used to launch CUDA's kernels
   */
  struct CUDAExecutionConfiguration {
    //! \brief number of threads per block
    int number_of_threads_per_block = 32;
  };

  /*!
   * \brief assign the evaluator to a function
   * \param[in] ctx: execution context
   * \param[in] c: CUDA execution configuration
   * \param[in] lhs: left hand side
   * \param[in] e: right hand side
   */
  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  [[nodiscard]] bool assign(AbstractErrorHandler&,
                            const CUDAExecutionConfiguration&,
                            FunctionType&,
                            const EvaluatorType)  //
      requires(
          ((LinearElementSpaceConcept<evaluator_space<EvaluatorType>>) ||
           (LinearQuadratureSpaceConcept<evaluator_space<EvaluatorType>>)) &&
          std::same_as<function_space<FunctionType>,
                       evaluator_space<EvaluatorType>>);

}  // end of namespace mgis::function

#include "MGIS/Function/CUDA/Algorithms.ixx"

#endif /* LIB_MGIS_FUNCTION_CUDA_ALGORITHMS_HXX */
