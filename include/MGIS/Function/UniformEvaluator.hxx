/*!
 * \file   MGIS/Function/UniformEvaluator.hxx
 * \brief File which defines the `UniformEvaluator` class, an evaluator which
 * returns the same value at each integration point of a space.
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
   * \brief class describing a uniform evaluator, i.e. an evaluator which
   * returns the same value at each integration point of a space
   *
   * The value is copied when the evaluator is built, so the object it is
   * built from may be destroyed afterwards.
   *
   * \tparam Space: type of the space on which the evaluator is defined
   * \tparam N: number of components of the returned values, which must be
   * greater than zero
   */
  template <SpaceConcept Space, std::size_t N>
  struct UniformEvaluator {
    /*!
     * \brief type of the value passed to construct the uniform evaluator,
     * which is a scalar value if `N` is equal to one, a constant
     * `std::span` of `N` values otherwise
     */
    using ValueType =
        std::conditional_t<N == 1, real, std::span<const real, N>>;
    /*!
     * \brief constructor
     * \param[in] s: space on which the evaluator is defined
     * \param[in] v: uniform value, which is copied in the evaluator
     */
    constexpr UniformEvaluator(const Space&, const ValueType&) MGIS_NOEXCEPT;
    //! \brief move constructor
    constexpr UniformEvaluator(UniformEvaluator&&) noexcept;
    //! \brief copy constructor
    constexpr UniformEvaluator(const UniformEvaluator&) MGIS_NOEXCEPT;
    //! \return the space on which the evaluator is defined
    [[nodiscard]] constexpr const Space& getSpace() const MGIS_NOEXCEPT;
    /*!
     * \brief call operator, which returns the uniform value associated with
     * an integration point
     * \param[in] i: index of the integration point, which is unused since
     * the value is the same everywhere
     * \return the uniform value
     * \note this overload is available when the space matches the
     * `ElementSpaceConcept` concept
     */
    [[nodiscard]] constexpr ValueType operator()(
        const element_index<Space>&) const MGIS_NOEXCEPT
      requires(ElementSpaceConcept<Space>);
    /*!
     * \brief call operator, which returns the uniform value associated with
     * an integration point
     * \param[in] wk: workspace associated with the element, which is unused
     * \param[in] i: index of the integration point, which is unused since
     * the value is the same everywhere
     * \return the uniform value
     * \note this overload is available when the space matches the
     * `ElementSpaceConcept` concept and exposes an element workspace
     */
    [[nodiscard]] constexpr ValueType operator()(
        const element_workspace<Space>&,
        const element_index<Space>&) const MGIS_NOEXCEPT
      requires(ElementSpaceConcept<Space> && hasElementWorkspace<Space>);
    /*!
     * \brief call operator, which returns the uniform value associated with
     * an integration point
     * \param[in] e: index of the cell, which is unused since the value is
     * the same everywhere
     * \param[in] i: index of the quadrature point, which is unused since
     * the value is the same everywhere
     * \return the uniform value
     * \note this overload is available when the space matches the
     * `QuadratureSpaceConcept` concept
     */
    [[nodiscard]] constexpr ValueType operator()(
        const cell_index<Space>&,
        const quadrature_point_index<Space>&) const MGIS_NOEXCEPT
      requires(QuadratureSpaceConcept<Space>);
    /*!
     * \brief call operator, which returns the uniform value associated with
     * an integration point
     * \param[in] wk: workspace associated with the cell, which is unused
     * \param[in] e: index of the cell, which is unused since the value is
     * the same everywhere
     * \param[in] i: index of the quadrature point, which is unused since
     * the value is the same everywhere
     * \return the uniform value
     * \note this overload is available when the space matches the
     * `QuadratureSpaceConcept` concept and exposes a cell workspace
     */
    [[nodiscard]] constexpr ValueType operator()(
        const cell_workspace<Space>&,
        const cell_index<Space>&,
        const quadrature_point_index<Space>&) const MGIS_NOEXCEPT
      requires(QuadratureSpaceConcept<Space> && hasCellWorkspace<Space>);
    //! \brief destructor
    constexpr ~UniformEvaluator() noexcept;

   protected:
    //! \brief a uniform evaluator has at least one component
    static_assert(N > 0);
    /*!
     * \brief builds the value stored by the evaluator from the scalar value
     * given to the constructor
     * \param[in] v: scalar value
     * \return the scalar value itself
     * \note this overload is used when `N` is equal to one
     */
    [[nodiscard]] static constexpr real buildValues(const real) MGIS_NOEXCEPT;
    /*!
     * \brief builds the array stored by the evaluator from the values given
     * to the constructor
     * \param[in] v: values, which must contain at least `N` values
     * \return an array holding the first `N` values of `v`
     * \note this overload is used when `N` is greater than one
     */
    [[nodiscard]] static constexpr std::array<real, N> buildValues(
        const std::span<const real>&) MGIS_NOEXCEPT;
    //! \brief underlying discretization space
    const Space space;
    //! \brief uniform value, a scalar if `N` is equal to one, an array
    //! otherwise
    const std::conditional_t<N == 1, real, std::array<real, N>> values;
  };  // end of struct UniformEvaluator

  /*!
   * \brief deduction guide: a `UniformEvaluator` built from a scalar value is
   * a scalar valued evaluator
   * \tparam SpaceType: type of the space on which the evaluator is defined
   */
  template <SpaceConcept SpaceType>
  UniformEvaluator(const SpaceType&, real) -> UniformEvaluator<SpaceType, 1>;

  /*!
   * \brief deduction guide: a `UniformEvaluator` built from a C array of `N`
   * values has `N` components
   * \tparam SpaceType: type of the space on which the evaluator is defined
   * \tparam N: number of components, which must be greater than one
   */
  template <SpaceConcept SpaceType, std::size_t N>
  UniformEvaluator(const SpaceType&, const real (&)[N])
      ->UniformEvaluator<SpaceType, N>
      /* requires(N > 1) */;  // most compilers (nvhpc, clang, intel) don't
                              // support require clause on class template
                              // deduction guide in 2025

  /*!
   * \brief deduction guide: a `UniformEvaluator` built from a `std::array` of
   * `N` values has `N` components
   * \tparam SpaceType: type of the space on which the evaluator is defined
   * \tparam N: number of components, which must be greater than one
   */
  template <SpaceConcept SpaceType, std::size_t N>
  UniformEvaluator(const SpaceType&, std::array<real, N>)
      ->UniformEvaluator<SpaceType, N>
      /* requires(N > 1) */;

  /*!
   * \brief deduction guide: a `UniformEvaluator` built from a `std::span` of
   * `N` values has `N` components
   * \tparam SpaceType: type of the space on which the evaluator is defined
   * \tparam N: number of components, which must be greater than one
   */
  template <SpaceConcept SpaceType, std::size_t N>
  UniformEvaluator(const SpaceType&, std::span<real, N>)
      ->UniformEvaluator<SpaceType, N>
      /* requires(N > 1) */;

  /*!
   * \brief checks that a uniform evaluator is consistent with the space on
   * which it is defined
   * \tparam Space: type of the space on which the evaluator is defined
   * \tparam N: number of components of the evaluator values
   * \param[in] eh: unused error handler
   * \param[in] e: uniform evaluator to be checked
   * \return always true: the values of a uniform evaluator are stored in the
   * evaluator itself and are always valid
   */
  template <SpaceConcept Space, std::size_t N>
  constexpr bool check(AbstractErrorHandler&,
                       const UniformEvaluator<Space, N>&) MGIS_NOEXCEPT;

  /*!
   * \brief return the space on which a uniform evaluator is defined
   * \tparam Space: type of the space on which the evaluator is defined
   * \tparam N: number of components of the evaluator values
   * \param[in] e: uniform evaluator
   * \return the space on which the evaluator is defined
   */
  template <SpaceConcept Space, std::size_t N>
  constexpr const Space& getSpace(const UniformEvaluator<Space, N>&)
      MGIS_NOEXCEPT;

  /*!
   * \brief return the number of components of a uniform evaluator
   * \tparam Space: type of the space on which the evaluator is defined
   * \tparam N: number of components of the evaluator values
   * \param[in] e: uniform evaluator whose number of components is queried
   * \return `N`, the number of components of the evaluator values
   */
  template <SpaceConcept Space, std::size_t N>
  constexpr mgis::size_type getNumberOfComponents(
      const UniformEvaluator<Space, N>&) MGIS_NOEXCEPT;

}  // end of namespace mgis::function

#include "MGIS/Function/UniformEvaluator.ixx"

#endif /* LIB_MGIS_FUNCTION_UNIFORMEVALUATOR_HXX */
