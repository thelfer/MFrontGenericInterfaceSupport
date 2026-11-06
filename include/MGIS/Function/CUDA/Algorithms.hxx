/*!
 * \file   MGIS/Function/CUDA/Algorithms.hxx
 * \brief  CUDA counterpart of `MGIS/Function/Algorithms.hxx`. It declares the
 * `assign` algorithm which evaluates an evaluator on the graphic card using the
 * CUDA programming model and the `CUDAExecutionConfiguration` type used to
 * configure the launch of the CUDA kernels.
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
   *
   * This is the CUDA counterpart of the `assign` algorithm declared in
   * `MGIS/Function/Algorithms.hxx`: the values are computed on the graphic card
   * by a CUDA kernel instead of being computed on the host. The evaluator and
   * the function must therefore be callable from the device, which means in
   * particular that they must be built on data accessible from the device, e.g.
   * `thrust::universal_vector`s, and that the lambda expressions they rely on
   * must be `MGIS_HOST_DEVICE` lambda expressions.
   *
   * The spaces of the function and of the evaluator and their number of
   * components are compared and the evaluator is checked before the kernel is
   * launched; the launch strategy itself is described in
   * `MGIS/Function/CUDA/Algorithms.ixx`.
   *
   * Only the function arguments which can be captured by a kernel are
   * supported: if the function is not a lightweight view, a view of it is
   * computed on the host and captured by the kernel.
   *
   * \tparam FunctionType: type of the function the values are assigned to,
   * which must be defined on the same space as the evaluator
   * \tparam EvaluatorType: type of the evaluator, which must satisfy the
   * `EvaluatorConcept` concept and be defined on a linear element space
   * \param[in] ctx: execution context, used to report the errors
   * \param[in] c: CUDA execution configuration
   * \param[in] f: left hand side, i.e. the function the values are assigned to
   * \param[in] e: right hand side, i.e. the evaluator computing the values
   * \return true if the assignment succeeded
   */
  template <typename FunctionType, EvaluatorConcept EvaluatorType>
  [[nodiscard]] bool assign(AbstractErrorHandler&,
                            const CUDAExecutionConfiguration&,
                            FunctionType&,
                            const EvaluatorType)  //
      requires(
          (LinearElementSpaceConcept<evaluator_space<EvaluatorType>>) &&
          std::same_as<function_space<FunctionType>,
                       evaluator_space<EvaluatorType>>);

}  // end of namespace mgis::function

#include "MGIS/Function/CUDA/Algorithms.ixx"

#endif /* LIB_MGIS_FUNCTION_CUDA_ALGORITHMS_HXX */
