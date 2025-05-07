/*!
 * \file   MGIS/Function/Algorithm.ixx
 * \brief
 * \author Thomas Helfer
 * \date   29/04/2025
 */

#ifndef LIB_MGIS_FUNCTION_ALGORITHMS_IXX
#define LIB_MGIS_FUNCTION_ALGORITHMS_IXX

#include <iterator>
#include <algorithm>

namespace mgis::function::algorithm {

  /*!
   * \brief copy N values
   * \tparam N: number of values to be copied
   * \tparam InputIterator: input iterator type
   * \tparam OutputIterator: output iterator type
   * \param[in] p: iterator to the beginning of the values to be copied
   * \param[in] pe: iterator past the end of the values to be copied
   * \param[in] po: iterator to the output values
   */
  template <size_type N, typename InputIterator, typename OutputIterator>
  void copy(const InputIterator p,
            const InputIterator pe,
            OutputIterator po) requires(N > 0) {
    if constexpr ((std::random_access_iterator<InputIterator>)&&  //
                  (std::random_access_iterator<OutputIterator>)) {
      if constexpr (N > 9) {
        std::copy(p, pe, po);
      } else if constexpr (N == 1) {
        *po = *p;
      } else if constexpr (N == 2) {
        po[0] = p[0];
        po[1] = p[1];
      } else if constexpr (N == 3) {
        po[0] = p[0];
        po[1] = p[1];
        po[2] = p[2];
      } else if constexpr (N == 4) {
        po[0] = p[0];
        po[1] = p[1];
        po[2] = p[2];
        po[3] = p[3];
      } else {
        copy<N - 1>(++p, pe, ++po);
      }
    } else {
      std::copy(p, pe, po);
    }
  }  // end of copy

}  // end of namespace mgis::function::algorithm

namespace mgis::function {

  template <size_type N, FunctionEvalutorConcept EvaluatorType>
  bool assign(Function& f, EvaluatorType e) requires(N > 0) {
    checkMatchingSpaces(f, e);
    raise_if(f.getNumberOfComponents() != N,
             "assign: invalid number of components for the left hand size");
    raise_if(e.getNumberOfComponents() != N,
             "assign: invalid number of components for the right hand size");
    //
    if (!e.check()) {
      return false;
    }
    e.allocateWorkspace();
    //
    const auto qspace = f.getSpace();
    const auto ne = qspace.size();
    for (size_type i = 0; i != ne; ++i) {
      if constexpr (N == 1) {
        auto& v = f.getValue(i);
        v = e(i);
      } else {
        auto lhs_values = f.template getValues<N>(i);
        const auto& rhs_values = e(i);
        algorithm::copy<N>(rhs_values.begin(), rhs_values.end(),
                           lhs_values.begin());
      }
    }
    return true;
  }  // end of assign

  template <FunctionEvalutorConcept EvaluatorType>
  bool assign(Function& f, EvaluatorType e) {
    raise_if(&f.getSpace() != &e.getSpace(),
             "assign: unmatched number of components for the left hand size "
             "and the right hand side");
    raise_if(f.getNumberOfComponents() != e.getNumberOfComponents(),
             "assign: unmatched number of components for the left hand size "
             "and the right hand side");
    //
    if (!e.check()) {
      return false;
    }
    e.allocateWorkspace();
    //
    const auto qspace = f.getSpace();
    const auto ne = qspace.size();
    if (f.isScalar()) {
      using result_type = std::invoke_result_t<EvaluatorType, size_type>;
      if constexpr (std::same_as<std::decay_t<result_type>, real>) {
        for (size_type i = 0; i != ne; ++i) {
          auto& lhs_value = f.getValue(i);
          lhs_value = e(i);
        }
      } else {
        for (size_type i = 0; i != ne; ++i) {
          auto lhs_value = f.getValue(i);
          lhs_value = *(e(i).begin());
        }
      }
    } else {
      for (size_type i = 0; i != ne; ++i) {
        auto lhs_values = f.getValues(i);
        const auto& rhs_values = e(i);
        std::copy(rhs_values.begin(), rhs_values.end(), lhs_values.begin());
      }
    }
    return true;
  }  // end of assign

}  // namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_ALGORITHMS_IXX */
