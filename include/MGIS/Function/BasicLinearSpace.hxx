/*!
 * \file   MGIS/Function/BasicLinearSpace.hxx
 * \brief
 * \author Thomas Helfer
 * \date   01/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_BASICLINEARSPACE_HXX
#define LIB_MGIS_FUNCTION_BASICLINEARSPACE_HXX

#include "MGIS/Config.hxx"
#include "MGIS/Function/SpaceConcept.hxx"

namespace mgis::function {

  //! \brief the simpliest implementation of a linear space
  struct MGIS_EXPORT BasicLinearSpace {
    static constexpr bool linear_element_indexing = true;
    using size_type = mgis::size_type;
    using element_index_type = mgis::size_type;
    /*!
     * \brief constructor
     * \param[in] s: size of the space
     */
    constexpr BasicLinearSpace(const size_type) noexcept;
    constexpr BasicLinearSpace(BasicLinearSpace&&) noexcept;
    constexpr BasicLinearSpace(const BasicLinearSpace&) noexcept;
    //
    constexpr size_type size() const noexcept;
    //! \brief destructor
    constexpr ~BasicLinearSpace() noexcept;

   private:
    //! \brief number of elements of the the space
    size_type nelts;
  };

}  // end of namespace mgis::function

#include "MGIS/Function/BasicLinearSpace.ixx"

#endif /* LIB_MGIS_FUNCTION_BASICLINEARSPACE_HXX */
