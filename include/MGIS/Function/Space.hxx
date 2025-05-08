/*!
 * \file   MGIS/Function/AbstractSpace.hxx
 * \brief
 * \author Thomas Helfer
 * \date   01/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_ABSTRACTSPACE_HXX
#define LIB_MGIS_FUNCTION_ABSTRACTSPACE_HXX

#include <concepts>
#include <type_traits>
#include "MGIS/Config.hxx"

namespace mgis::function {

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
      using type = typename Space::element_index_type;
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
      using type = typename Space::ElementWorkspace;
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
      using type = typename Space::cell_index_type;
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
      using type = typename Space::quadrature_point_index_type;
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
      using type = typename Space::CellWorkspace;
    };

  }  // namespace internals

  //! \brief a concept describing a space
  template <typename SpaceType>
  concept SpaceConcept = requires(const SpaceType& s) {
    std::integral<typename SpaceType::size_type>;
    { s.size() } -> internals::same_decay_type<typename SpaceType::size_type>;
  };

  template <typename SpaceType>
  concept ElementSpaceConcept = requires(const SpaceType& s) {
    SpaceConcept<SpaceType>;
    typename SpaceType::element_index_type;
  };

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
    typename SpaceType::ElementWorkspace;
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
   * are stored for 0 to size() - 1
   */
  template <typename SpaceType>
  concept LinearElementSpaceConcept = requires(const SpaceType& s) {
    ElementSpaceConcept<SpaceType>;
    SpaceType::linear_element_indexing;
    std::integral<typename SpaceType::element_index_type>;
  };

  /*!
   * \brief a concept describing a quadrature space where
   * values can be addressed by given the element index
   * and the value of the integration points
   */
  template <typename SpaceType>
  concept QuadratureSpaceConcept = requires(const SpaceType& s) {
    typename SpaceType::cell_index_type;
    std::integral<typename SpaceType::quadrature_point_index_type>;
    { s.getNumberOfCells() } -> std::same_as<typename SpaceType::size_type>;
    {
      s.getNumberOfQuadraturePoints(
          std::declval<typename SpaceType::cell_index_type>())
      } -> internals::same_decay_type<
          typename SpaceType::quadrature_point_index_type>;
  };

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
  concept LinearQuadratureSpaceConcept = QuadratureSpaceConcept<SpaceType> &&
      requires(const SpaceType& s) {
    SpaceType::linear_cell_indexing;
    std::integral<typename SpaceType::cell_index_type>;
    {
      s.getQuadraturePointOffset(
          std::declval<typename SpaceType::cell_index_type>(),
          std::declval<typename SpaceType::quadrature_point_index_type>())
      } -> std::integral;
  };

  /*!
   * \brief  a simple meta function stating if SpaceType
   * does not expose an CellWorkspace structure
   */
  template <SpaceConcept SpaceType>
  inline constexpr bool hasCellWorkspace = requires {
    typename SpaceType::CellWorkspace;
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

#endif /* LIB_MGIS_FUNCTION_ABSTRACTSPACE_HXX */
