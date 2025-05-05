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
    using size_type = Invalid;
    using ElementWorkspace = Invalid;
    using element_index_type = Invalid;
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
  concept LinearSpaceConcept = requires(const SpaceType& s) {
    std::integral<typename SpaceTraits<std::decay_t<SpaceType>>::size_type>;
    {
      s.size()
      }
      -> std::same_as<typename SpaceTraits<std::decay_t<SpaceType>>::size_type>;
  };

  /*!
   * \brief a concept describing a quadrature space where
   * values can be addressed by given the element index
   * and the value of the integration points
   */
  template <typename SpaceType>
  concept QuadratureSpaceConcept = requires(const SpaceType& s) {
    !std::same_as<
        typename SpaceTraits<std::decay_t<SpaceType>>::ElementWorkspace,
        typename SpaceTraitsBase::Invalid>;
    !std::same_as<
        typename SpaceTraits<std::decay_t<SpaceType>>::element_index_type,
        typename SpaceTraitsBase::Invalid>;
    std::integral<typename SpaceTraits<
        std::decay_t<SpaceType>>::quadrature_point_index_type>;
    {
      s.getNumberOfElements()
      } -> std::same_as<
          typename SpaceTraits<std::decay_t<SpaceType>>::element_index_type>;
    {
      s.getNumberOfQuadraturePoints(
          std::declval<typename SpaceTraits<
              std::decay_t<SpaceType>>::element_index_type>())
      } -> std::same_as<typename SpaceTraits<
          std::decay_t<SpaceType>>::quadrature_point_index_type>;
    {
      s.getElementWorkspace(
          std::declval<typename SpaceTraits<
              std::decay_t<SpaceType>>::element_index_type>())
      } -> std::same_as<
          typename SpaceTraits<std::decay_t<SpaceType>>::ElementWorkspace>;
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
      LinearSpaceConcept<SpaceType> || QuadratureSpaceConcept<SpaceType>;

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_ABSTRACTSPACE_HXX */
