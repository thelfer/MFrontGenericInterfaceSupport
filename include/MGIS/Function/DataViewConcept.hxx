/*!
 * \file   MGIS/Function/DataViewConcept.hxx
 * \brief  This file declares `DataViewConcept`
 * \author Thomas Helfer
 * \date   02/10/2026
 */

#ifndef LIB_MGIS_FUNCTION_DATAVIEWCONCEPT_HXX
#define LIB_MGIS_FUNCTION_DATAVIEWCONCEPT_HXX

#include <span>
#include <type_traits>

namespace mgis::function::internals {

  /*!
   * \brief this class must be specialised
   * for class holding views to data, such as
   * std::span or tfel::math::View
   */
  template <typename T>
  struct IsDataView : std::false_type {};

  //! \brief partial specialisation for std::span
  template <class T, std::size_t Extent>
  struct IsDataView<::std::span<T, Extent>> : std::true_type{};

}  // end of namespace mgis::function::internals

namespace mgis::function {

  //! \brief a concept used to distinguish views from data
  template <typename T>
  concept DataViewConcept = ::mgis::function::internals::IsDataView<T>::type;

}  // end of namespace mgis::function

#endif /* LIB_MGIS_FUNCTION_DATAVIEWCONCEPT_HXX */
