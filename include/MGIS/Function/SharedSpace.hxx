/*!
 * \file   MGIS/Function/SharedSpace.hxx
 * \brief Definition of `mgis::function::SharedSpace`, a class which
 * encapsulates a space in a `std::shared_ptr` so that this space can be shared
 * by several functions or evaluators, and of the functions which query the
 * encapsulated space.
 * \author Thomas Helfer
 * \date   02/06/2025
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MFEM_MGIS_SHAREDSPACE_HXX
#define LIB_MFEM_MGIS_SHAREDSPACE_HXX

#include <memory>
#include "MGIS/Raise.hxx"
#include "MGIS/Contract.hxx"
#include "MGIS/Function/SpaceConcept.hxx"

namespace mgis::function::internals {

  /*!
   * \brief helper class used to build the traits of a `SharedSpace`: this
   * primary template adds no `ElementWorkspace` type, which is the case when
   * the underlying space does not expose such a type
   *
   * \tparam HasElementWorkspace: states if the underlying space exposes an
   * `ElementWorkspace` type
   * \tparam SpaceType: type of the underlying space
   */
  template <bool, ElementSpaceConcept>
  struct SharedElementSpaceAddElementWorkspace {};

  /*!
   * \brief specialization used when the underlying space exposes an
   * `ElementWorkspace` type: this type is forwarded to the traits of the shared
   * space
   *
   * \tparam SpaceType: type of the underlying space
   */
  template <ElementSpaceConcept SpaceType>
  struct SharedElementSpaceAddElementWorkspace<true, SpaceType> {
    //! \brief the element workspace type of the underlying space
    using ElementWorkspace = typename SpaceTraits<SpaceType>::ElementWorkspace;
  };

  /*!
   * \brief traits of a `SharedSpace` whose underlying space does not match the
   * `ElementSpaceConcept` concept, which are thus empty
   *
   * \tparam SpaceType: type of the underlying space
   */
  template <typename>
  struct SharedElementSpaceTraits {};

  /*!
   * \brief specialization used when the underlying space matches the
   * `ElementSpaceConcept` concept: the index type and the element workspace
   * type of the underlying space are forwarded to the traits of the shared
   * space
   *
   * \tparam SpaceType: type of the underlying space
   */
  template <ElementSpaceConcept SpaceType>
  struct SharedElementSpaceTraits<SpaceType>
      : SharedElementSpaceAddElementWorkspace<hasElementWorkspace<SpaceType>,
                                              SpaceType> {
    //! \brief the element index type of the underlying space
    using element_index_type =
        typename SpaceTraits<SpaceType>::element_index_type;
  };

  /*!
   * \brief helper class used to build the traits of a `SharedSpace`: this
   * primary template adds nothing, as the underlying space does not match the
   * `LinearElementSpaceConcept` concept
   *
   * \tparam SpaceType: type of the underlying space
   */
  template <typename>
  struct SharedLinearElementSpaceTraits {};

  /*!
   * \brief specialization used when the underlying space matches the
   * `LinearElementSpaceConcept` concept: the elements of the shared space are
   * addressed exactly as those of the underlying space
   *
   * \tparam SpaceType: type of the underlying space
   */
  template <LinearElementSpaceConcept SpaceType>
  struct SharedLinearElementSpaceTraits<SpaceType> {
    /*!
     * \brief the elements of the space are addressed by an index varying from 0
     * to the size of the space minus one
     */
    static constexpr bool linear_element_indexing = true;
  };

  /*!
   * \brief helper class used to build the traits of a `SharedSpace`: this
   * primary template adds no `CellWorkspace` type, which is the case when the
   * underlying space does not expose such a type
   *
   * \tparam HasCellWorkspace: states if the underlying space exposes a
   * `CellWorkspace` type
   * \tparam SpaceType: type of the underlying space
   */
  template <bool, QuadratureSpaceConcept>
  struct SharedQuadratureSpaceAddCellWorkspace {};

  /*!
   * \brief specialization used when the underlying space exposes a
   * `CellWorkspace` type: this type is forwarded to the traits of the shared
   * space
   *
   * \tparam SpaceType: type of the underlying space
   */
  template <QuadratureSpaceConcept SpaceType>
  struct SharedQuadratureSpaceAddCellWorkspace<true, SpaceType> {
    //! \brief the cell workspace type of the underlying space
    using CellWorkspace = typename SpaceTraits<SpaceType>::CellWorkspace;
  };

  /*!
   * \brief traits of a `SharedSpace` whose underlying space does not match the
   * `QuadratureSpaceConcept` concept, which are thus empty
   *
   * \tparam SpaceType: type of the underlying space
   */
  template <typename>
  struct SharedQuadratureSpaceTraits {};

  /*!
   * \brief specialization used when the underlying space matches the
   * `QuadratureSpaceConcept` concept: the index types of the underlying space
   * are forwarded to the traits of the shared space
   *
   * \tparam SpaceType: type of the underlying space
   */
  template <QuadratureSpaceConcept SpaceType>
  struct SharedQuadratureSpaceTraits<SpaceType>
      : SharedQuadratureSpaceAddCellWorkspace<hasCellWorkspace<SpaceType>,
                                              SpaceType> {
    //! \brief the cell index type of the underlying space
    using cell_index_type = typename SpaceTraits<SpaceType>::cell_index_type;
    //! \brief the quadrature point index type of the underlying space
    using quadrature_point_index_type =
        typename SpaceTraits<SpaceType>::quadrature_point_index_type;
  };

  /*!
   * \brief helper class used to build the traits of a `SharedSpace`: this
   * primary template adds nothing, as the underlying space does not match the
   * `LinearQuadratureSpaceConcept` concept
   *
   * \tparam SpaceType: type of the underlying space
   */
  template <typename>
  struct SharedLinearQuadratureSpaceTraits {};

  /*!
   * \brief specialization used when the underlying space matches the
   * `LinearQuadratureSpaceConcept` concept: the cells of the shared space are
   * addressed exactly as those of the underlying space
   *
   * \tparam SpaceType: type of the underlying space
   */
  template <LinearQuadratureSpaceConcept SpaceType>
  struct SharedLinearQuadratureSpaceTraits<SpaceType> {
    /*!
     * \brief the cells of the space are addressed by an index varying from 0 to
     * the number of cells minus one
     */
    static constexpr bool linear_cell_indexing = true;
  };

}  // end of namespace mgis::function::internals

namespace mgis::function {

  /*!
   * \brief this class provides a way to encapsulate a space into a sharable
   * object
   *
   * Such an object behaves exactly as the underlying space (it satisfies the
   * same concepts) except for the `areEquivalent` function, which only states
   * that two shared spaces are equivalent when they point to the same object.
   *
   * \tparam SpaceType: type of the space encapsulated in the sharable object
   */
  template <SpaceConcept SpaceType>
  struct SharedSpace : private PreconditionsChecker<SharedSpace<SpaceType>> {
    /*!
     * \brief checks the preconditions of the constructor, i.e. that the given
     * pointer is not null
     * \param[in,out] eh: error handler used to report the errors
     * \param[in] s: pointer to the space to encapsulate
     * \return false if the given pointer is null, true otherwise
     * \note the "invalid space" error message is registered in the error
     * handler when the given pointer is null
     */
    [[nodiscard]] static bool checkPreconditions(
        AbstractErrorHandler& eh, const std::shared_ptr<const SpaceType>& s) {
      if (s.get() == nullptr) {
        return eh.registerErrorMessage("invalid space");
      }
      return true;
    }  // end of checkPreconditions

    /*!
     * \brief constructor, which checks the preconditions before building the
     * object
     * \param[in] s: pointer to the space to encapsulate
     * \note a contract violation is reported if the given pointer is null
     */
    SharedSpace(const std::shared_ptr<const SpaceType>& s)
        : SharedSpace(preconditions_check, s) {}

    /*!
     * \brief constructor
     * \tparam doPreconditionChecks: if true, the preconditions are checked,
     * otherwise they are skipped
     * \param[in] pcheck: tag stating if the preconditions are checked
     * (`preconditions_check`) or skipped (`no_precondition_check`)
     * \param[in] s: pointer to the space to encapsulate
     */
    template <bool doPreconditionChecks>
    SharedSpace(const PreconditionsCheck<doPreconditionChecks>& pcheck,
                const std::shared_ptr<const SpaceType>& s)
        : PreconditionsChecker<SharedSpace>(pcheck, s),
          space(s) {}  // end of SharedSpace

    //! \brief move constructor
    SharedSpace(SharedSpace&&) = default;
    //! \brief copy constructor
    SharedSpace(const SharedSpace&) = default;

    /*!
     * \return a constant reference to the encapsulated space
     */
    const SpaceType& operator*() const MGIS_NOEXCEPT { return *(this->space); }
    /*!
     * \return a constant pointer to the encapsulated space
     */
    const SpaceType* get() const MGIS_NOEXCEPT { return this->space.get(); }

   private:
    //! \brief pointer to the encapsulated space
    std::shared_ptr<const SpaceType> space;
  };

  // template deduction guide lines

  /*!
   * \brief deduction guide for a pointer to a constant space
   * \tparam SpaceType: type of the space pointed to
   */
  template <SpaceConcept SpaceType>
  SharedSpace(const std::shared_ptr<const SpaceType>&)
      -> SharedSpace<SpaceType>;

  /*!
   * \brief deduction guide for a pointer to a space which may be modified
   * \tparam SpaceType: type of the space pointed to
   */
  template <SpaceConcept SpaceType>
  SharedSpace(const std::shared_ptr<SpaceType>&) -> SharedSpace<SpaceType>;

  /*!
   * \brief specialization of the `SpaceTraits` class for a shared space: the
   * properties of the underlying space which make sense for a shared space are
   * forwarded by the helpers of the `mgis::function::internals` namespace
   *
   * \tparam SpaceType: type of the space encapsulated in the sharable object
   */
  template <SpaceConcept SpaceType>
  struct SpaceTraits<SharedSpace<SpaceType>>
      : internals::SharedElementSpaceTraits<SpaceType>,
        internals::SharedLinearElementSpaceTraits<SpaceType>,
        internals::SharedQuadratureSpaceTraits<SpaceType>,
        internals::SharedLinearQuadratureSpaceTraits<SpaceType> {
    //! \brief type used to store the size of the space
    using size_type = typename SpaceTraits<SpaceType>::size_type;
  };

  /*!
   * \brief checks if two shared spaces are equivalent
   * \tparam SpaceType: type of the encapsulated spaces
   * \param[in] s: first shared space
   * \param[in] s2: second shared space
   * \return true if the two shared spaces point to the same object
   * \note contrary to the `areEquivalent` function of the underlying space,
   * this function does not compare the properties of the spaces but the
   * identity of the objects which own them
   */
  template <SpaceConcept SpaceType>
  [[nodiscard]] bool areEquivalent(const SharedSpace<SpaceType>& s,
                                   const SharedSpace<SpaceType>& s2)
      MGIS_NOEXCEPT {
    return s.get() == s2.get();
  }

  /*!
   * \brief returns the size of the encapsulated space
   * \tparam SpaceType: type of the encapsulated space
   * \param[in] s: shared space
   * \return the number of elements of the encapsulated space
   */
  template <SpaceConcept SpaceType>
  [[nodiscard]] auto getSpaceSize(const SharedSpace<SpaceType>& s) {
    return getSpaceSize(*s);
  }

  /*!
   * \brief returns the number of cells of the encapsulated space
   * \tparam SpaceType: type of the encapsulated space
   * \param[in] s: shared space
   * \return the number of cells of the encapsulated space
   */
  template <QuadratureSpaceConcept SpaceType>
  [[nodiscard]] auto getNumberOfCells(const SharedSpace<SpaceType>& s) {
    return getNumberOfCells(*s);
  }

  /*!
   * \brief returns the number of quadrature points of a given cell of the
   * encapsulated space
   * \tparam SpaceType: type of the encapsulated space
   * \param[in] s: shared space
   * \param[in] e: cell index
   * \return the number of quadrature points of the cell of index `e`
   */
  template <QuadratureSpaceConcept SpaceType>
  [[nodiscard]] auto getNumberOfQuadraturePoints(
      const SharedSpace<SpaceType>& s, const cell_index<SpaceType>& e) {
    return getNumberOfQuadraturePoints(*s, e);
  }

  /*!
   * \brief returns the offset of a given quadrature point of the encapsulated
   * space
   * \tparam SpaceType: type of the encapsulated space
   * \param[in] s: shared space
   * \param[in] e: cell index
   * \param[in] i: quadrature point index of the cell of index `e`
   * \return the offset of the considered quadrature point
   */
  template <LinearQuadratureSpaceConcept SpaceType>
  [[nodiscard]] auto getQuadraturePointOffset(
      const SharedSpace<SpaceType>& s,
      const cell_index<SpaceType>& e,
      const quadrature_point_index<SpaceType>& i) {
    return getQuadraturePointOffset(*s, e, i);
  }

  /*!
   * \brief returns the workspace associated with a given element of the
   * encapsulated space
   * \tparam SpaceType: type of the encapsulated space
   * \param[in] s: shared space
   * \param[in] i: element index
   * \return the workspace of the element of index `i`
   * \note this function is only defined if the encapsulated space exposes an
   * `ElementWorkspace` type
   */
  template <ElementSpaceConcept SpaceType>
  [[nodiscard]] decltype(auto) getElementWorkspace(
      const SharedSpace<SpaceType>& s, const element_index<SpaceType>& i)  //
      requires(hasElementWorkspace<SpaceType>) {
    return getElementWorkspace(*s, i);
  }

  /*!
   * \brief returns the workspace associated with a given cell of the
   * encapsulated space
   * \tparam SpaceType: type of the encapsulated space
   * \param[in] s: shared space
   * \param[in] e: cell index
   * \return the workspace of the cell of index `e`
   * \note this function is only defined if the encapsulated space exposes a
   * `CellWorkspace` type
   */
  template <QuadratureSpaceConcept SpaceType>
  [[nodiscard]] decltype(auto) getCellWorkspace(
      const SharedSpace<SpaceType>& s, const cell_index<SpaceType>& e)  //
      requires(hasCellWorkspace<SpaceType>) {
    return getCellWorkspace(*s, e);
  }

}  // end of namespace mgis::function

#endif /* LIB_MFEM_MGIS_SHAREDSPACE_HXX */
