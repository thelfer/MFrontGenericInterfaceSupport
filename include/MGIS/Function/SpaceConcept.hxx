/*!
 * \file   MGIS/Function/SpaceConcept.hxx
 * \brief File which defines the concept a space must satisfy and the helpers
 * used to query the index types and the workspaces of a space.
 * \author Thomas Helfer
 * \date   01/05/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_FUNCTION_SPACECONCEPT_HXX
#define LIB_MGIS_FUNCTION_SPACECONCEPT_HXX

#include <concepts>
#include <type_traits>
#include "MGIS/Config.hxx"

namespace mgis::function {

  /*!
   * \brief a quadrature traits class defining the main properties of a space
   *
   * For a type to match the `SpaceConcept`,  `SpaceTraits<SpaceType>` must
   * expose a type alias named `size_type`.
   *
   * For a type to match the `ElementSpaceConcept`,  `SpaceTraits<SpaceType>`
   * must expose a type alias named `element_index_type`.
   *
   * For a type to match the `LinearElementSpaceConcept`, the
   * `element_index_type` type alias of `SpaceTraits<SpaceType>` must match the
   * `std::integral` concept and `SpaceTraits<SpaceType>` must expose a
   * `static constexpr`
   * boolean data member named `linear_element_indexing` that must be equal to
   * `true`.
   *
   * For a type to match the `QuadratureSpaceConcept`,  `SpaceTraits<SpaceType>`
   * must expose a type alias named `cell_index_type` and another type alias
   * named `quadrature_point_index_type` which must match the `std::integral`
   * concept.
   *
   * For a type to match the `LinearQuadratureSpaceConcept`, the
   * `cell_index_type` type alias of `SpaceTraits<SpaceType>` must match the
   * `std::integral`
   * concept and `SpaceTraits<SpaceType>` must expose a `static constexpr`
   * boolean data member named `linear_cell_indexing` that must be equal to
   * `true`.
   *
   * \tparam SpaceType: type of the space considered
   */
  template <typename SpaceType>
  struct SpaceTraits {};

  namespace internals {

    /*!
     * \brief concept satisfied if the given types are the same once
     * `std::decay_t` is applied
     *
     * \tparam T: first type
     * \tparam U: second type
     */
    template <typename T, typename U>
    concept same_decay_type = std::same_as<std::decay_t<T>, std::decay_t<U>>;

    /*!
     * \brief helper metafunction returning an undefined type, used as fallback
     * when a space does not expose the queried type
     */
    struct UndefinedTypeSelector {
      //! \brief placeholder standing for an undefined type
      struct Undefined;
      //! \brief the undefined type
      using type = Undefined;
    };

    /*!
     * \brief metafunction used to return the `element_index_type` type alias
     * of `SpaceTraits<Space>` if it is defined,
     * `UndefinedTypeSelector::Undefined` otherwise
     *
     * \tparam hasElementIndexType: boolean stating if the
     * `element_index_type` type alias of `SpaceTraits<Space>` is defined
     * \tparam Space: type of the space considered
     */
    template <bool hasElementIndexType, typename Space>
    struct ElementIndexTypeSelector : UndefinedTypeSelector {};

    /*!
     * \brief specialization used when the space exposes an `element_index_type`
     * type alias
     *
     * \tparam Space: type of the space considered
     */
    template <typename Space>
    struct ElementIndexTypeSelector<true, Space> {
      //! \brief the element index type exposed by the space
      using type = typename SpaceTraits<Space>::element_index_type;
    };

    /*!
     * \brief metafunction used to return the `ElementWorkspace` type of
     * `SpaceTraits<Space>` if it is defined,
     * `UndefinedTypeSelector::Undefined` otherwise
     *
     * \tparam hasElementWorkspace: boolean stating if the
     * `ElementWorkspace` type of `SpaceTraits<Space>` is defined
     * \tparam Space: type of the space considered
     */
    template <bool hasElementWorkspace, typename Space>
    struct ElementWorkspaceSelector : UndefinedTypeSelector {};

    /*!
     * \brief specialization used when the space exposes an `ElementWorkspace`
     * type
     *
     * \tparam Space: type of the space considered
     */
    template <typename Space>
    struct ElementWorkspaceSelector<true, Space> {
      //! \brief the element workspace type exposed by the space
      using type = typename SpaceTraits<Space>::ElementWorkspace;
    };

    /*!
     * \brief metafunction used to return the `cell_index_type` type alias of
     * `SpaceTraits<Space>` if it is defined,
     * `UndefinedTypeSelector::Undefined` otherwise
     *
     * \tparam hasCellIndexType: boolean stating if the `cell_index_type` type
     * alias of `SpaceTraits<Space>` is defined
     * \tparam Space: type of the space considered
     */
    template <bool hasCellIndexType, typename Space>
    struct CellIndexTypeSelector : UndefinedTypeSelector {};

    /*!
     * \brief specialization used when the space exposes a `cell_index_type`
     * type alias
     *
     * \tparam Space: type of the space considered
     */
    template <typename Space>
    struct CellIndexTypeSelector<true, Space> {
      //! \brief the cell index type exposed by the space
      using type = typename SpaceTraits<Space>::cell_index_type;
    };

    /*!
     * \brief metafunction used to return the
     * `quadrature_point_index_type` type alias of `SpaceTraits<Space>` if it
     * is defined, `UndefinedTypeSelector::Undefined` otherwise
     *
     * \tparam hasQuadraturePointIndexType: boolean stating if the
     * `quadrature_point_index_type` type alias of `SpaceTraits<Space>` is
     * defined
     * \tparam Space: type of the space considered
     */
    template <bool hasQuadraturePointIndexType, typename Space>
    struct QuadraturePointIndexTypeSelector : UndefinedTypeSelector {};

    /*!
     * \brief specialization used when the space exposes a
     * `quadrature_point_index_type` type alias
     *
     * \tparam Space: type of the space considered
     */
    template <typename Space>
    struct QuadraturePointIndexTypeSelector<true, Space> {
      //! \brief the quadrature point index type exposed by the space
      using type = typename SpaceTraits<Space>::quadrature_point_index_type;
    };

    /*!
     * \brief metafunction used to return the `CellWorkspace` type of
     * `SpaceTraits<Space>` if it is defined,
     * `UndefinedTypeSelector::Undefined` otherwise
     *
     * \tparam hasCellWorkspace: boolean stating if the `CellWorkspace` type of
     * `SpaceTraits<Space>` is defined
     * \tparam Space: type of the space considered
     */
    template <bool hasCellWorkspace, typename Space>
    struct CellWorkspaceSelector : UndefinedTypeSelector {};

    /*!
     * \brief specialization used when the space exposes a `CellWorkspace` type
     *
     * \tparam Space: type of the space considered
     */
    template <typename Space>
    struct CellWorkspaceSelector<true, Space> {
      //! \brief the cell workspace type exposed by the space
      using type = typename SpaceTraits<Space>::CellWorkspace;
    };

  }  // namespace internals

  /*!
   * \brief a concept describing a space
   *
   * A space is a type which:
   * - exposes a `size_type` type alias through `SpaceTraits`, which must match
   *   the `std::integral` concept,
   * - can be moved constructed,
   * - provides its size through the `getSpaceSize` function,
   * - states if two spaces are equivalent through the `areEquivalent` function.
   *
   * \tparam SpaceType: type of the space considered
   */
  template <typename SpaceType>
  concept SpaceConcept =
      ((std::integral<typename SpaceTraits<SpaceType>::size_type>)&&  //
       (std::move_constructible<SpaceType>)&&                         //
       (requires(const SpaceType& s) {
         {
           getSpaceSize(s)
           } -> internals::same_decay_type<
               typename SpaceTraits<SpaceType>::size_type>;
       }) &&
       (requires(const SpaceType& s, const SpaceType& s2) {
         { areEquivalent(s, s2) } -> std::same_as<bool>;
       }));

  /*!
   * \brief a concept describing a space where values can be addressed by an
   * element index
   *
   * \tparam SpaceType: type of the space considered
   */
  template <typename SpaceType>
  concept ElementSpaceConcept =
      ((SpaceConcept<SpaceType>)&&  //
       requires { typename SpaceTraits<SpaceType>::element_index_type; });

  /*!
   * \brief a simple alias giving the type of the index of an element of a space
   *
   * \tparam SpaceType: type of the space considered
   * \note this alias leads to an undefined type if `SpaceType` does not match
   * the `ElementSpaceConcept` concept
   */
  template <SpaceConcept SpaceType>
  using element_index = typename internals::
      ElementIndexTypeSelector<ElementSpaceConcept<SpaceType>, SpaceType>::type;

  /*!
   * \brief a simple meta function stating if `SpaceType` exposes an
   * `ElementWorkspace` type
   *
   * \tparam SpaceType: type of the space considered
   */
  template <SpaceConcept SpaceType>
  inline constexpr bool hasElementWorkspace = requires {
    typename SpaceTraits<SpaceType>::ElementWorkspace;
  };

  /*!
   * \brief a simple alias giving the type of the workspace associated with an
   * element of a space
   *
   * \tparam SpaceType: type of the space considered
   * \note this alias leads to an undefined type if `SpaceType` does not expose
   * an `ElementWorkspace` type
   */
  template <SpaceConcept SpaceType>
  using element_workspace = typename internals::
      ElementWorkspaceSelector<hasElementWorkspace<SpaceType>, SpaceType>::type;

  /*!
   * \brief a concept describing a space where all the elements are stored from
   * 0 to `getSpaceSize(s) - 1`
   *
   * \tparam SpaceType: type of the space considered
   */
  template <typename SpaceType>
  concept LinearElementSpaceConcept =
      ((ElementSpaceConcept<SpaceType>) &&                                 //
       (requires { SpaceTraits<SpaceType>::linear_element_indexing; }) &&  //
       (SpaceTraits<SpaceType>::linear_element_indexing) &&                //
       (std::integral<
           typename SpaceTraits<SpaceType>::element_index_type>) &&  //
       (std::same_as<typename SpaceTraits<SpaceType>::size_type,
                     typename SpaceTraits<SpaceType>::element_index_type>));

  /*!
   * \brief a concept describing a quadrature space where values can be
   * addressed given the cell index and the index of the integration points
   *
   * \tparam SpaceType: type of the space considered
   */
  template <typename SpaceType>
  concept QuadratureSpaceConcept =
      ((requires(const SpaceType& s) {
         typename SpaceTraits<SpaceType>::cell_index_type;
         typename SpaceTraits<SpaceType>::quadrature_point_index_type;
         {
           getNumberOfCells(s)
           } -> std::same_as<typename SpaceTraits<SpaceType>::size_type>;
         {
           getNumberOfQuadraturePoints(
               s,
               std::declval<typename SpaceTraits<SpaceType>::cell_index_type>())
           } -> internals::same_decay_type<
               typename SpaceTraits<SpaceType>::quadrature_point_index_type>;
       }) &&  //
       (std::integral<
           typename SpaceTraits<SpaceType>::quadrature_point_index_type>));

  /*!
   * \brief a simple alias giving the type of the index of a cell of a space
   *
   * \tparam SpaceType: type of the space considered
   * \note this alias leads to an undefined type if `SpaceType` does not match
   * the `QuadratureSpaceConcept` concept
   */
  template <SpaceConcept SpaceType>
  using cell_index = typename internals::
      CellIndexTypeSelector<QuadratureSpaceConcept<SpaceType>, SpaceType>::type;

  /*!
   * \brief a simple alias giving the type of the index of a quadrature point of
   * a space
   *
   * \tparam SpaceType: type of the space considered
   * \note this alias leads to an undefined type if `SpaceType` does not match
   * the `QuadratureSpaceConcept` concept
   */
  template <SpaceConcept SpaceType>
  using quadrature_point_index =
      typename internals::QuadraturePointIndexTypeSelector<
          QuadratureSpaceConcept<SpaceType>,
          SpaceType>::type;

  /*!
   * \brief a concept describing a quadrature space where all the cells are
   * stored from 0 to `getNumberOfCells(s) - 1` and the offset of a quadrature
   * point can be computed from the index of its cell and its own index
   *
   * \tparam SpaceType: type of the space considered
   */
  template <typename SpaceType>
  concept LinearQuadratureSpaceConcept =
      ((QuadratureSpaceConcept<SpaceType>)&&  //
       (requires { SpaceTraits<SpaceType>::linear_cell_indexing; }) &&
       (SpaceTraits<SpaceType>::linear_cell_indexing) &&
       (std::integral<typename SpaceTraits<SpaceType>::cell_index_type>)&&  //
       (requires(const SpaceType& s) {
         {
           getQuadraturePointOffset(
               s,
               std::declval<typename SpaceTraits<SpaceType>::cell_index_type>(),
               std::declval<typename SpaceTraits<
                   SpaceType>::quadrature_point_index_type>())
           } -> internals::same_decay_type<
               typename SpaceTraits<SpaceType>::size_type>;
       }));

  /*!
   * \brief a simple meta function stating if `SpaceType` exposes a
   * `CellWorkspace` type
   *
   * \tparam SpaceType: type of the space considered
   */
  template <SpaceConcept SpaceType>
  inline constexpr bool hasCellWorkspace = requires {
    typename SpaceTraits<SpaceType>::CellWorkspace;
  };

  /*!
   * \brief a simple alias giving the type of the workspace associated with a
   * cell of a space
   *
   * \tparam SpaceType: type of the space considered
   * \note this alias leads to an undefined type if `SpaceType` does not expose
   * a `CellWorkspace` type
   */
  template <SpaceConcept SpaceType>
  using cell_workspace =
      typename internals::CellWorkspaceSelector<hasCellWorkspace<SpaceType>,
                                                SpaceType>::type;

  /*!
   * \brief a concept describing a space where values can be addressed either by
   * an element index or by a cell index and a quadrature point index
   *
   * \tparam SpaceType: type of the space considered
   */
  template <typename SpaceType>
  concept FunctionalSpaceConcept =
      ElementSpaceConcept<SpaceType> || QuadratureSpaceConcept<SpaceType>;

  /*!
   * \brief a concept describing a space which matches both the
   * `FunctionalSpaceConcept` concept and either the `LinearElementSpaceConcept`
   * concept or the `LinearQuadratureSpaceConcept` concept
   *
   * \tparam Space: type of the space considered
   */
  template <typename Space>
  concept LinearFunctionalSpaceConcept = ((FunctionalSpaceConcept<Space>)&&(
      LinearElementSpaceConcept<Space> || LinearQuadratureSpaceConcept<Space>));

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_SPACECONCEPT_HXX */
