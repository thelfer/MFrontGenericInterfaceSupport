/*!
 * \file   MGIS/Function/UniformEvaluator.ixx
 * \brief
 * \author Thomas Helfer
 * \date   06/10/2025
 */

#ifndef LIB_MGIS_FUNCTION_UNIFORMEVALUATOR_IXX
#define LIB_MGIS_FUNCTION_UNIFORMEVALUATOR_IXX

#include <utility>

namespace mgis::function {

  template <SpaceConcept Space, std::size_t N>
  constexpr real UniformEvaluator<Space, N>::buildValues(
      const real v) noexcept {
    return v;
  }  // end of buildValues

  template <SpaceConcept Space, std::size_t N>
  constexpr std::array<real, N> UniformEvaluator<Space, N>::buildValues(
      const std::span<const real>& v) noexcept {
    return [&v]<std::size_t... Is>(std::index_sequence<Is...>) {
      return std::array<real, N>{{v[Is]...}};
    }
    (std::make_index_sequence<N>());
  }  // end of buildValues

  template <SpaceConcept Space, std::size_t N>
  constexpr UniformEvaluator<Space, N>::UniformEvaluator(
      const Space& s, const ValueType& v) noexcept
      : space(s), values(buildValues(v)) {}  // end of UniformEvaluator

  template <SpaceConcept Space, std::size_t N>
  constexpr UniformEvaluator<Space, N>::UniformEvaluator(
      UniformEvaluator&&) noexcept = default;

  template <SpaceConcept Space, std::size_t N>
  constexpr UniformEvaluator<Space, N>::UniformEvaluator(
      const UniformEvaluator&) noexcept = default;

  template <SpaceConcept Space, std::size_t N>
  constexpr const Space& UniformEvaluator<Space, N>::getSpace() const noexcept {
    return this->space;
  }  // end of space

  template <SpaceConcept Space, std::size_t N>
  constexpr typename UniformEvaluator<Space, N>::ValueType
  UniformEvaluator<Space, N>::operator()(const element_index<Space>&)
      const noexcept requires(ElementSpaceConcept<Space>) {
    return this->values;
  }  // end of operator()

  template <SpaceConcept Space, std::size_t N>
  constexpr typename UniformEvaluator<Space, N>::ValueType
  UniformEvaluator<Space, N>::operator()(
      const element_workspace<Space>&,
      const element_index<Space>&) const noexcept
      requires(ElementSpaceConcept<Space>&& hasElementWorkspace<Space>) {
    return this->values;
  }  // end of operator()

  template <SpaceConcept Space, std::size_t N>
  constexpr typename UniformEvaluator<Space, N>::ValueType
  UniformEvaluator<Space, N>::operator()(const cell_index<Space>&,
                                         const quadrature_point_index<Space>&)
      const noexcept requires(QuadratureSpaceConcept<Space>) {
    return this->values;
  }  // end of operator()

  template <SpaceConcept Space, std::size_t N>
  constexpr typename UniformEvaluator<Space, N>::ValueType
  UniformEvaluator<Space, N>::operator()(
      const cell_workspace<Space>&,
      const cell_index<Space>&,
      const quadrature_point_index<Space>&) const noexcept
      requires(QuadratureSpaceConcept<Space>&& hasCellWorkspace<Space>) {
    return this->values;
  }  // end of operator()

  template <SpaceConcept Space, std::size_t N>
  constexpr UniformEvaluator<Space, N>::~UniformEvaluator() noexcept = default;

  template <SpaceConcept Space, std::size_t N>
  constexpr bool check(AbstractErrorHandler&,
                       const UniformEvaluator<Space, N>&) noexcept {
    return true;
  }  // end of check

  template <SpaceConcept Space, std::size_t N>
  constexpr const Space& getSpace(
      const UniformEvaluator<Space, N>& e) noexcept {
    return e.getSpace();
  }  // end of getSpace

  template <SpaceConcept Space, std::size_t N>
  constexpr mgis::size_type getNumberOfComponents(
      const UniformEvaluator<Space, N>&) noexcept {
    return static_cast<mgis::size_type>(N);
  }  // end of getSpace

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_UNIFORMEVALUATOR_IXX */
