/*!
 * \file   MGIS/Function/TFEL/DataView.hxx
 * \brief  This file specializes the `IsDataView` class for TFEL' types
 * \author Thomas Helfer
 * \date   02/10/2026
 */

#ifndef LIB_MGIS_FUNCTION_TFEL_DATAVIEW_HXX
#define LIB_MGIS_FUNCTION_TFEL_DATAVIEW_HXX

#include "TFEL/Math/Array/View.hxx"
#include "TFEL/Math/Array/CoalescedView.hxx"
#include "MGIS/Function/DataViewConcept.hxx"

namespace mgis::function::internals {

  //! \brief partial specialisation for tfel::math::View
  template <typename MappedType, typename IndexingPolicyType>
  struct IsDataView<::tfel::math::View<MappedType, IndexingPolicyType>>
    : std::true_type{};

  //! \brief partial specialisation for tfel::math::CoalescedView
  template <
      ::tfel::math::MappableMathObjectUsingCoalescedViewConcept MappedType,
      typename IndexingPolicyType>
  struct IsDataView<::tfel::math::CoalescedView<MappedType, IndexingPolicyType>>
    : std::true_type{};

}  // end of namespace mgis::function::internals

#endif /* LIB_MGIS_FUNCTION_TFEL_DATAVIEW_HXX */
