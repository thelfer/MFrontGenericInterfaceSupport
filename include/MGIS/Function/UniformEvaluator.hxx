/*!
 * \file   MGIS/Function/UniformEvaluator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   06/10/2025
 */

#ifndef LIB_MGIS_FUNCTION_UNIFORMEVALUATOR_HXX
#define LIB_MGIS_FUNCTION_UNIFORMEVALUATOR_HXX

#include <span>
#include <array>
#include <type_traits>
#include "MGIS/Config.hxx"
#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/EvaluatorConcept.hxx"

namespace mgis::function {

  /*!
   * \brief class describing a uniform evaluator
   */
  template <SpaceConcept Space, std::size_t N>
  struct UniformEvaluator {
    //! \brief type passed to construct the uniform value
    using ValueType =
        std::conditional_t<N == 1, real, std::span<const real, N>>;
    /*!
     * \brief constructor
     * \param[in] s: space
     * \param[in] v: value
     */
    constexpr UniformEvaluator(const Space&, const ValueType&) noexcept;
    //! \brief move constructor
    constexpr UniformEvaluator(UniformEvaluator&&) noexcept;
    //! \brief copy constructor
    constexpr UniformEvaluator(const UniformEvaluator&) noexcept;
    //! \brief return the underlying space
    constexpr const Space& getSpace() const noexcept;
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    constexpr ValueType operator()(const element_index<Space>&) const noexcept
        requires(ElementSpaceConcept<Space>);
    /*!
     * \brief call operator
     * \param[in] i: integration point index
     */
    constexpr ValueType operator()(const element_workspace<Space>&,
                                   const element_index<Space>&) const noexcept
        requires(ElementSpaceConcept<Space>&& hasElementWorkspace<Space>);
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    constexpr ValueType operator()(const cell_index<Space>&,
                                   const quadrature_point_index<Space>&)
        const noexcept requires(QuadratureSpaceConcept<Space>);
    /*!
     * \brief call operator
     * \param[in] e: cell index
     * \param[in] i: integration point index
     */
    constexpr ValueType operator()(
        const cell_workspace<Space>&,
        const cell_index<Space>&,
        const quadrature_point_index<Space>&) const noexcept
        requires(QuadratureSpaceConcept<Space>&& hasCellWorkspace<Space>);
    //! \brief destructor
    constexpr ~UniformEvaluator() noexcept;

   protected:
    //
    static_assert(N > 0);
    //
    static constexpr real buildValues(const real) noexcept;
    //
    static constexpr std::array<real, N> buildValues(
        const std::span<const real>&) noexcept;
    //! \brief underlying discretization space
    const Space space;
    //! \brief values
    const std::conditional_t<N == 1, real, std::array<real, N>> values;
  };  // end of struct UniformEvaluator

  // class template deduction guide

  template <SpaceConcept SpaceType>
  UniformEvaluator(const SpaceType&, real) -> UniformEvaluator<SpaceType, 1>;

  template <SpaceConcept SpaceType, std::size_t N>
  UniformEvaluator(const SpaceType&, const real (&)[N])
      ->UniformEvaluator<SpaceType, N>
  requires(N > 1);

  template <SpaceConcept SpaceType, std::size_t N>
  UniformEvaluator(const SpaceType&, std::array<real, N>)
      ->UniformEvaluator<SpaceType, N>
  requires(N > 1);

  template <SpaceConcept SpaceType, std::size_t N>
  UniformEvaluator(const SpaceType&, std::span<real, N>)
      ->UniformEvaluator<SpaceType, N>
  requires(N > 1);

  template <SpaceConcept Space, std::size_t N>
  constexpr bool check(AbstractErrorHandler&,
                       const UniformEvaluator<Space, N>&) noexcept;

  template <SpaceConcept Space, std::size_t N>
  constexpr const Space& getSpace(const UniformEvaluator<Space, N>&) noexcept;

  template <SpaceConcept Space, std::size_t N>
  constexpr mgis::size_type getNumberOfComponents(
      const UniformEvaluator<Space, N>&) noexcept;

}  // end of namespace mgis::function

#include "MGIS/Function/UniformEvaluator.ixx"

#endif /* LIB_MGIS_FUNCTION_UNIFORMEVALUATOR_HXX */
