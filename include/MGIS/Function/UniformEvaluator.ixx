/*!
 * \file   MGIS/Function/UniformEvaluator.ixx
 * \brief Inline implementation of the classes declared in
 * `MGIS/Function/UniformEvaluator.hxx`, which includes this file at its end.
 * \author Thomas Helfer
 * \date   06/10/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_UNIFORMEVALUATOR_IXX
#define LIB_MGIS_FUNCTION_UNIFORMEVALUATOR_IXX

#include <utility>

namespace mgis::function {

  template <SpaceConcept Space, std::size_t N>
  constexpr real UniformEvaluator<Space, N>::buildValues(
      const real v) MGIS_NOEXCEPT {
    return v;
  }  // end of buildValues

  template <SpaceConcept Space, std::size_t N>
  constexpr std::array<real, N> UniformEvaluator<Space, N>::buildValues(
      const std::span<const real>& v) MGIS_NOEXCEPT {
    // the `std::index_sequence` is expanded so that the first `N` values of
    // the span are copied in the array stored by the evaluator
    return [&v]<std::size_t... Is>(std::index_sequence<Is...>) {
      return std::array<real, N>{{v[Is]...}};
    }
    (std::make_index_sequence<N>());
  }  // end of buildValues

  template <SpaceConcept Space, std::size_t N>
  constexpr UniformEvaluator<Space, N>::UniformEvaluator(
      const Space& s, const ValueType& v) MGIS_NOEXCEPT
      : space(s), values(buildValues(v)) {}  // end of UniformEvaluator

  template <SpaceConcept Space, std::size_t N>
  constexpr UniformEvaluator<Space, N>::UniformEvaluator(
      UniformEvaluator&&) noexcept = default;

  template <SpaceConcept Space, std::size_t N>
  constexpr UniformEvaluator<Space, N>::UniformEvaluator(
      const UniformEvaluator&) MGIS_NOEXCEPT = default;

  template <SpaceConcept Space, std::size_t N>
  constexpr const Space& UniformEvaluator<Space, N>::getSpace() const MGIS_NOEXCEPT {
    return this->space;
  }  // end of space

  template <SpaceConcept Space, std::size_t N>
  constexpr typename UniformEvaluator<Space, N>::ValueType
  UniformEvaluator<Space, N>::operator()(const element_index<Space>&)
      const MGIS_NOEXCEPT requires(ElementSpaceConcept<Space>) {
    return this->values;
  }  // end of operator()

  template <SpaceConcept Space, std::size_t N>
  constexpr typename UniformEvaluator<Space, N>::ValueType
  UniformEvaluator<Space, N>::operator()(
      const element_workspace<Space>&,
      const element_index<Space>&) const MGIS_NOEXCEPT
      requires(ElementSpaceConcept<Space>&& hasElementWorkspace<Space>) {
    return this->values;
  }  // end of operator()

  template <SpaceConcept Space, std::size_t N>
  constexpr typename UniformEvaluator<Space, N>::ValueType
  UniformEvaluator<Space, N>::operator()(const cell_index<Space>&,
                                         const quadrature_point_index<Space>&)
      const MGIS_NOEXCEPT requires(QuadratureSpaceConcept<Space>) {
    return this->values;
  }  // end of operator()

  template <SpaceConcept Space, std::size_t N>
  constexpr typename UniformEvaluator<Space, N>::ValueType
  UniformEvaluator<Space, N>::operator()(
      const cell_workspace<Space>&,
      const cell_index<Space>&,
      const quadrature_point_index<Space>&) const MGIS_NOEXCEPT
      requires(QuadratureSpaceConcept<Space>&& hasCellWorkspace<Space>) {
    return this->values;
  }  // end of operator()

  template <SpaceConcept Space, std::size_t N>
  constexpr UniformEvaluator<Space, N>::~UniformEvaluator() noexcept = default;

  template <SpaceConcept Space, std::size_t N>
  constexpr bool check(AbstractErrorHandler&,
                       const UniformEvaluator<Space, N>&) MGIS_NOEXCEPT {
    return true;
  }  // end of check

  template <SpaceConcept Space, std::size_t N>
  constexpr const Space& getSpace(
      const UniformEvaluator<Space, N>& e) MGIS_NOEXCEPT {
    return e.getSpace();
  }  // end of getSpace

  template <SpaceConcept Space, std::size_t N>
  constexpr mgis::size_type getNumberOfComponents(
      const UniformEvaluator<Space, N>&) MGIS_NOEXCEPT {
    return static_cast<mgis::size_type>(N);
  }  // end of getSpace

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_UNIFORMEVALUATOR_IXX */
