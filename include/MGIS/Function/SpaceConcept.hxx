/*!
 * \file   MGIS/Function/SpaceConcept.hxx
 * \brief
 * \author Thomas Helfer
 * \date   01/05/2025
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
   * For a type to match the `LinearElementSpaceConcept`,
   * `SpaceTraits<SpaceType>::element_index_type` must match the `std::integral`
   * concept and `SpaceTraits<SpaceType>` must expose a `static constexpr`
   * boolean data member named `linear_element_indexing` that must be equal to
   * `true`.
   *
   * For a type to match the `QuadratureSpaceConcept`,  `SpaceTraits<SpaceType>`
   * must expose a type alias named `cell_index_type` and another type alias
   * named `quadrature_point_index_type` which must match the `std::integral`
   * concept.
   *
   * For a type to match the `LinearQuadratureSpaceConcept`,
   `SpaceTraits<SpaceType>::cell_index_type` must match the `std::integral`
   * concept and `SpaceTraits<SpaceType>` must expose a `static constexpr`
   * boolean data member named `linear_cell_indexing` that must be equal to
   * `true`.
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
     * \brief a metafunction returning an undefined type
     */
    struct UndefinedTypeSelector {
      //! \brief a simple place holder
      struct Undefined;
      //! \brief result
      using type = Undefined;
    };

    /*!
     * \brief a meta function used to return SpaceType::element_index_type if it
     * is defined, UndefinedTypeSelector::Undefined otherwise
     *
     * \tparam hasElementIndexType: boolean stating if
     * SpaceType::ElementIndexType exists
     * \tparam SpaceType: type of the space
     * considered
     */
    template <bool hasElementIndexType, typename Space>
    struct ElementIndexTypeSelector : UndefinedTypeSelector {};

    template <typename Space>
    struct ElementIndexTypeSelector<true, Space> {
      //! \brief the type exposed by Space
      using type = typename SpaceTraits<Space>::element_index_type;
    };

    /*!
     * \brief a meta function used to return SpaceType::ElementWorkspace if it
     * is defined, UndefinedTypeSelector::Undefined otherwise
     *
     * \tparam hasElementWorkspace: boolean stating if
     * SpaceType::ElementWorkspace exists
     * \tparam SpaceType: type of the space
     * considered
     */
    template <bool hasElementWorkspace, typename Space>
    struct ElementWorkspaceSelector : UndefinedTypeSelector {};

    template <typename Space>
    struct ElementWorkspaceSelector<true, Space> {
      //! \brief the type exposed by Space
      using type = typename SpaceTraits<Space>::ElementWorkspace;
    };

    /*!
     * \brief a meta function used to return SpaceType::cell_index_type if it
     * is defined, UndefinedTypeSelector::Undefined otherwise
     *
     * \tparam hasCellIndexType: boolean stating if
     * SpaceType::CellIndexType exists
     * \tparam SpaceType: type of the space
     * considered
     */
    template <bool hasCellIndexType, typename Space>
    struct CellIndexTypeSelector : UndefinedTypeSelector {};

    template <typename Space>
    struct CellIndexTypeSelector<true, Space> {
      //! \brief the type exposed by Space
      using type = typename SpaceTraits<Space>::cell_index_type;
    };

    /*!
     * \brief a meta function used to return
     * SpaceType::quadrature_point_index_type if it is defined,
     * UndefinedTypeSelector::Undefined otherwise
     *
     * \tparam hasQuadraturePointIndexType: boolean stating if
     * SpaceType::QuadraturePointIndexType exists
     * \tparam SpaceType: type of the space
     * considered
     */
    template <bool hasQuadraturePointIndexType, typename Space>
    struct QuadraturePointIndexTypeSelector : UndefinedTypeSelector {};

    template <typename Space>
    struct QuadraturePointIndexTypeSelector<true, Space> {
      //! \brief the type exposed by Space
      using type = typename SpaceTraits<Space>::quadrature_point_index_type;
    };

    /*!
     * \brief a meta function used to return SpaceType::CellWorkspace if it
     * is defined, UndefinedTypeSelector::Undefined otherwise
     *
     * \tparam hasCellWorkspace: boolean stating if
     * SpaceType::CellWorkspace exists
     * \tparam SpaceType: type of the space
     * considered
     */
    template <bool hasCellWorkspace, typename Space>
    struct CellWorkspaceSelector : UndefinedTypeSelector {};

    template <typename Space>
    struct CellWorkspaceSelector<true, Space> {
      //! \brief the type exposed by Space
      using type = typename SpaceTraits<Space>::CellWorkspace;
    };

  }  // namespace internals

  //! \brief a concept describing a space
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

  template <typename SpaceType>
  concept ElementSpaceConcept =
      ((SpaceConcept<SpaceType>)&&  //
       requires { typename SpaceTraits<SpaceType>::element_index_type; });

  //! \brief a simple alias
  template <SpaceConcept SpaceType>
  using element_index = typename internals::
      ElementIndexTypeSelector<ElementSpaceConcept<SpaceType>, SpaceType>::type;

  /*!
   * \brief  a simple meta function stating if SpaceType
   * does not expose an ElementWorkspace structure
   */
  template <ElementSpaceConcept SpaceType>
  inline constexpr bool hasElementWorkspace = requires {
    typename SpaceTraits<SpaceType>::ElementWorkspace;
  };

  /*!
   * \brief a simple alias
   * \note this alias leads to an undefined type if SpaceType
   * does not expose an ElementWorkspace structure
   */
  template <SpaceConcept SpaceType>
  using element_workspace = typename internals::
      ElementWorkspaceSelector<hasElementWorkspace<SpaceType>, SpaceType>::type;

  /*!
   * \brief a concept describing a space where all the elements
   * are stored for 0 to getSpaceSize(s) - 1
   */
  template <typename SpaceType>
  concept LinearElementSpaceConcept =
      ((ElementSpaceConcept<SpaceType>)&&                                  //
       (requires { SpaceTraits<SpaceType>::linear_element_indexing; }) &&  //
       (std::integral<
           typename SpaceTraits<SpaceType>::element_index_type>)&&  //
       (std::same_as<typename SpaceTraits<SpaceType>::size_type,
                     typename SpaceTraits<SpaceType>::element_index_type>));

  /*!
   * \brief a concept describing a quadrature space where
   * values can be addressed by given the element index
   * and the value of the integration points
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

  //! \brief a simple alias
  template <SpaceConcept SpaceType>
  using cell_index = typename internals::
      CellIndexTypeSelector<QuadratureSpaceConcept<SpaceType>, SpaceType>::type;

  //! \brief a simple alias
  template <SpaceConcept SpaceType>
  using quadrature_point_index =
      typename internals::QuadraturePointIndexTypeSelector<
          QuadratureSpaceConcept<SpaceType>,
          SpaceType>::type;

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
   * \brief  a simple meta function stating if SpaceType
   * does not expose an CellWorkspace structure
   */
  template <SpaceConcept SpaceType>
  inline constexpr bool hasCellWorkspace = requires {
    typename SpaceTraits<SpaceType>::CellWorkspace;
  };

  /*!
   * \brief a simple alias
   * \note this alias leads to an undefined type if SpaceType
   * does not expose an CellWorkspace structure
   */
  template <SpaceConcept SpaceType>
  using cell_workspace =
      typename internals::CellWorkspaceSelector<hasCellWorkspace<SpaceType>,
                                                SpaceType>::type;

  template <typename SpaceType>
  concept FunctionalSpaceConcept =
      ElementSpaceConcept<SpaceType> || QuadratureSpaceConcept<SpaceType>;

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_SPACECONCEPT_HXX */
