/*!
 * \file   MGIS/Function/BasicLinearSpace.hxx
 * \brief Definition of `mgis::function::BasicLinearSpace`, the simplest
 * implementation of a space where the elements are addressed by an index
 * varying from 0 to the size of the space minus one.
 * \author Thomas Helfer
 * \date   01/05/2025
 */

#ifndef LIB_MGIS_FUNCTION_BASICLINEARSPACE_HXX
#define LIB_MGIS_FUNCTION_BASICLINEARSPACE_HXX

#include "MGIS/Config.hxx"
#include "MGIS/Function/SpaceConcept.hxx"

namespace mgis::function {

  /*!
   * \brief the simplest implementation of a linear space
   *
   * The elements of such a space are addressed by an index varying from 0 to
   * `getSpaceSize(s) - 1`.
   */
  struct BasicLinearSpace {
    /*!
     * \brief constructor
     * \param[in] s: size of the space
     */
    constexpr BasicLinearSpace(const size_type s) MGIS_NOEXCEPT;
    //! \brief move constructor
    constexpr BasicLinearSpace(BasicLinearSpace&&) noexcept;
    //! \brief copy constructor
    constexpr BasicLinearSpace(const BasicLinearSpace&) MGIS_NOEXCEPT;
    /*!
     * \return the size of the space, i.e. the number of elements it contains
     */
    [[nodiscard]] constexpr size_type size() const MGIS_NOEXCEPT;
    //! \brief destructor
    constexpr ~BasicLinearSpace() noexcept;

   private:
    //! \brief number of elements of the space
    size_type nelts;
  };

  //! \brief traits defining the properties of a `BasicLinearSpace`
  template <>
  struct SpaceTraits<BasicLinearSpace> {
    /*!
     * \brief the elements of the space are addressed by an index varying from 0
     * to the size of the space minus one
     */
    static constexpr bool linear_element_indexing = true;
    //! \brief type used to store the size of the space
    using size_type = mgis::size_type;
    //! \brief type of the index of an element of the space
    using element_index_type = mgis::size_type;
  };

  /*!
   * \brief returns the size of the space
   * \param[in] s: space
   * \return the number of elements the space contains
   */
  constexpr size_type getSpaceSize(const BasicLinearSpace& s) MGIS_NOEXCEPT;

  /*!
   * \brief checks if two spaces are equivalent, i.e. if they have the same size
   * \param[in] s1: first space
   * \param[in] s2: second space
   * \return true if the two spaces are equivalent
   */
  [[nodiscard]] constexpr bool areEquivalent(
      const BasicLinearSpace& s1, const BasicLinearSpace& s2) MGIS_NOEXCEPT;

}  // end of namespace mgis::function

#include "MGIS/Function/BasicLinearSpace.ixx"

#endif /* LIB_MGIS_FUNCTION_BASICLINEARSPACE_HXX */
