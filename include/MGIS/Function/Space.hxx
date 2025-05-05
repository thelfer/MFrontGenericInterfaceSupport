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

  /*!
   * \brief a base class for all traits class related to spaces.
   */
  struct SpaceTraitsBase {
    struct Invalid;
    // type that shall be returned by the size method
    using size_type = Invalid;
    //
    static constexpr bool linear_element_indexing = false;
    using element_index_type = Invalid;
    //
    static constexpr bool linear_cell_indexing = false;
    using CellWorkspace = Invalid;
    using cell_index_type = Invalid;
    using quadrature_point_index_type = Invalid;
  };

  /*!
   * \brief a traits class describing a space
   * \tparam SpaceType: type of the space
   */
  template <typename SpaceType>
  struct SpaceTraits : SpaceTraitsBase {};

  /*!
   * \brief a concept describing a space where all the elements
   * are stored for 0 to size() - 1
   */
  template <typename SpaceType>
  concept SpaceConcept = requires(const SpaceType& s) {
    std::integral<typename SpaceTraits<std::decay_t<SpaceType>>::size_type>;
    {
      s.size()
      }
      -> std::same_as<typename SpaceTraits<std::decay_t<SpaceType>>::size_type>;
  };

  /*!
   * \brief a concept describing a space where all the elements
   * are stored for 0 to size() - 1
   */
  template <typename SpaceType>
  concept LinearSpaceConcept = requires(const SpaceType& s) {
    SpaceConcept<SpaceType>;
    SpaceTraits<SpaceType>::linear_element_indexing;
    SpaceTraits<std::decay_t<SpaceType>>::linear_element_indexing
        ? std::integral<
              typename SpaceTraits<std::decay_t<SpaceType>>::element_index_type>
        : true;
  };

  /*!
   * \brief a concept describing a quadrature space where
   * values can be addressed by given the element index
   * and the value of the integration points
   */
  template <typename SpaceType>
  concept QuadratureSpaceConcept = requires(const SpaceType& s) {
    !std::same_as<
        typename SpaceTraits<std::decay_t<SpaceType>>::CellWorkspace,
        typename SpaceTraitsBase::Invalid>;
    !std::same_as<
        typename SpaceTraits<std::decay_t<SpaceType>>::cell_index_type,
        typename SpaceTraitsBase::Invalid>;
    SpaceTraits<std::decay_t<SpaceType>>::linear_cell_indexing
        ? std::integral<
              typename SpaceTraits<std::decay_t<SpaceType>>::cell_index_type>
        : true;
    std::integral<typename SpaceTraits<
        std::decay_t<SpaceType>>::quadrature_point_index_type>;
    {
      s.getNumberOfCells()
      }
      -> std::same_as<typename SpaceTraits<std::decay_t<SpaceType>>::size_type>;
    {
      s.getNumberOfQuadraturePoints(
          std::declval<typename SpaceTraits<
              std::decay_t<SpaceType>>::cell_index_type>())
      } -> std::same_as<typename SpaceTraits<
          std::decay_t<SpaceType>>::quadrature_point_index_type>;
    {
      s.getCellWorkspace(std::declval<typename SpaceTraits<
                                std::decay_t<SpaceType>>::cell_index_type>())
      } -> std::same_as<
          typename SpaceTraits<std::decay_t<SpaceType>>::CellWorkspace>;
  };

  template <typename SpaceType>
  concept LinearQuadratureSpaceConcept = LinearSpaceConcept<SpaceType> &&
      QuadratureSpaceConcept<SpaceType> &&  //
      requires(const SpaceType& s) {
    {
      s.getQuadraturePointOffset(
          std::declval<typename SpaceTraits<
              std::decay_t<SpaceType>>::element_index_type>(),
          std::declval<typename SpaceTraits<
              std::decay_t<SpaceType>>::quadrature_point_index_type>())
      } -> std::integral;
  };

  template <typename SpaceType>
  concept FunctionalSpaceConcept =
      SpaceConcept<SpaceType> || QuadratureSpaceConcept<SpaceType>;

  template <typename SpaceType>
  concept LinearFunctionalSpaceConcept =
      LinearSpaceConcept<SpaceType> || LinearQuadratureSpaceConcept<SpaceType>;

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_ABSTRACTSPACE_HXX */
