/*!
 * \file   MGIS/Function/DataViewConcept.hxx
 * \brief File which defines the concept used to distinguish the types holding
 * a view to some values, such as `std::span`, from the types owning those
 * values.
 * \author Thomas Helfer
 * \date   02/10/2026
 */

#ifndef LIB_MGIS_FUNCTION_DATAVIEWCONCEPT_HXX
#define LIB_MGIS_FUNCTION_DATAVIEWCONCEPT_HXX

#include <span>
#include <type_traits>

namespace mgis::function::internals {

  /*!
   * \brief traits stating if the given type is a view to some values, which
   * must be specialized for the classes holding views to data, such as
   * `std::span` or `tfel::math::View`
   *
   * The primary template evaluates to `std::false_type`.
   *
   * \tparam T: type to be tested
   */
  template <typename T>
  struct IsDataView : std::false_type {};

  /*!
   * \brief specialization for reference types: a reference does not own the
   * value it refers to, so it is always considered as a view
   *
   * \tparam T: referenced type
   */
  template <typename T>
  struct IsDataView<T&> : std::true_type {};

  /*!
   * \brief specialization for `std::span`
   *
   * \tparam T: type of the values pointed to by the span
   * \tparam Extent: number of values pointed to by the span, which may be
   * `std::dynamic_extent`
   */
  template <class T, std::size_t Extent>
  struct IsDataView<::std::span<T, Extent>> : std::true_type {};

}  // end of namespace mgis::function::internals

namespace mgis::function {

  /*!
   * \brief concept which is satisfied by the types holding a view to some
   * values owned by another object
   *
   * This concept is used to know if the result of a function or an evaluator
   * can be remapped without copying the values it refers to.
   *
   * \tparam T: type to be tested
   */
  template <typename T>
  concept DataViewConcept = ::mgis::function::internals::IsDataView<T>::value;

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_DATAVIEWCONCEPT_HXX */
