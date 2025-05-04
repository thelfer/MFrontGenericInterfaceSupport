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
   * \brief a traits class describing a space
   * \tparam SpaceType: type of the space
   */
  template <typename SpaceType>
  struct SpaceTraits;

  /*!
   * \brief a concept describing a space where all the elements
   * are stored for 0 to size() - 1
   */
  template <typename SpaceType>
  concept LinearSpaceConcept = requires(const SpaceType& s) {
    std::integral<typename SpaceTraits<SpaceType>::size_type>;
    { s.size() } -> std::same_as<typename SpaceTraits<SpaceType>::size_type>;
  };

  /*!
   * \brief a concept describing a quadrature space where
   * values can be addressed by given the element index
   * and the value of the integration points
   */
  template <typename SpaceType>
  concept QuadratureSpaceConcept = requires(const SpaceType& s) {
    std::integral<typename SpaceTraits<SpaceType>::quadrature_point_index_type>;
    {
      s.getNumberOfElements()
      } -> std::same_as<typename SpaceTraits<SpaceType>::element_index_type>;
    {
      s.getNumberOfQuadraturePoints(
          std::declval<typename SpaceTraits<SpaceType>::element_index_type>())
      } -> std::same_as<typename SpaceTraits<SpaceType>::quadrature_point_index_type>;
  };

  template <typename SpaceType>
  concept FunctionalSpaceConcept =
      LinearSpaceConcept<SpaceType> || QuadratureSpaceConcept<SpaceType>;

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_ABSTRACTSPACE_HXX */
