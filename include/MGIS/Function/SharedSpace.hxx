/*!
 * \file   MGIS/Function/SharedSpace.hxx
 * \brief
 * \author Thomas Helfer
 * \date   02/06/2025
 */

#ifndef LIB_MFEM_MGIS_SHAREDSPACE_HXX
#define LIB_MFEM_MGIS_SHAREDSPACE_HXX

#include <memory>
#include "MGIS/Raise.hxx"
#include "MGIS/Contract.hxx"
#include "MGIS/Function/SpaceConcept.hxx"

namespace mgis::function::internals {

  template <bool, ElementSpaceConcept>
  struct SharedElementSpaceAddElementWorkspace {};

  template <ElementSpaceConcept SpaceType>
  struct SharedElementSpaceAddElementWorkspace<true, SpaceType> {
    using ElementWorkspace = typename SpaceTraits<SpaceType>::ElementWorkspace;
  };

  template <typename>
  struct SharedElementSpaceTraits {};

  template <ElementSpaceConcept SpaceType>
  struct SharedElementSpaceTraits<SpaceType>
      : SharedElementSpaceAddElementWorkspace<hasElementWorkspace<SpaceType>,
                                              SpaceType> {
    using element_index_type =
        typename SpaceTraits<SpaceType>::element_index_type;
  };

  template <typename>
  struct SharedLinearElementSpaceTraits {};

  template <LinearElementSpaceConcept SpaceType>
  struct SharedLinearElementSpaceTraits<SpaceType> {
    static constexpr bool linear_element_indexing = true;
  };

  template <bool, QuadratureSpaceConcept>
  struct SharedQuadratureSpaceAddCellWorkspace {};

  template <QuadratureSpaceConcept SpaceType>
  struct SharedQuadratureSpaceAddCellWorkspace<true, SpaceType> {
    using CellWorkspace = typename SpaceTraits<SpaceType>::CellWorkspace;
  };

  template <typename>
  struct SharedQuadratureSpaceTraits {};

  template <QuadratureSpaceConcept SpaceType>
  struct SharedQuadratureSpaceTraits<SpaceType>
      : SharedQuadratureSpaceAddCellWorkspace<hasCellWorkspace<SpaceType>,
                                              SpaceType> {
    using cell_index_type = typename SpaceTraits<SpaceType>::cell_index_type;
    using quadrature_point_index_type =
        typename SpaceTraits<SpaceType>::quadrature_point_index_type;
  };

  template <typename>
  struct SharedLinearQuadratureSpaceTraits {};

  template <LinearQuadratureSpaceConcept SpaceType>
  struct SharedLinearQuadratureSpaceTraits<SpaceType> {
    static constexpr bool linear_cell_indexing = true;
  };

}  // end of namespace mgis::function::internals

namespace mgis::function {

  template <SpaceConcept SpaceType>
  struct SharedSpace : private PreconditionsChecker<SharedSpace<SpaceType>> {
    /*!
     * \param[in] eh: error handler
     */
    [[nodiscard]] static bool checkPreconditions(
        AbstractErrorHandler& eh, const std::shared_ptr<const SpaceType>& s) {
      if (s.get() == nullptr) {
        return eh.registerErrorMessage("invalid space");
      }
      return true;
    }  // end of checkPreconditions

    SharedSpace(const std::shared_ptr<const SpaceType>& s)
        : SharedSpace(preconditions_check, s){}

    template <bool doPreconditionChecks>
    SharedSpace(const PreconditionsCheck<doPreconditionChecks>& pcheck,
                const std::shared_ptr<const SpaceType>& s)
        : PreconditionsChecker<SharedSpace>(pcheck, s),
          space(s) {}  // end of SharedSpace

    SharedSpace(SharedSpace&&) = default;
    SharedSpace(const SharedSpace&) = default;

    const SpaceType& operator*() const noexcept { return *(this->space); }
    const SpaceType* get() const noexcept { return this->space.get(); }

   private:
    std::shared_ptr<const SpaceType> space;
  };

  template <SpaceConcept SpaceType>
  SharedSpace(const std::shared_ptr<const SpaceType>&)
      -> SharedSpace<SpaceType>;

  template <SpaceConcept SpaceType>
  SharedSpace(const std::shared_ptr<SpaceType>&) -> SharedSpace<SpaceType>;

  template <SpaceConcept SpaceType>
  struct SpaceTraits<SharedSpace<SpaceType>>
      : internals::SharedElementSpaceTraits<SpaceType>,
        internals::SharedLinearElementSpaceTraits<SpaceType>,
        internals::SharedQuadratureSpaceTraits<SpaceType>,
        internals::SharedLinearQuadratureSpaceTraits<SpaceType> {
    using size_type = typename SpaceTraits<SpaceType>::size_type;
  };

  template <SpaceConcept SpaceType>
  bool areEquivalent(const SharedSpace<SpaceType>& s,
                     const SharedSpace<SpaceType>& s2) noexcept {
    return s.get() == s2.get();
  }

  template <SpaceConcept SpaceType>
  auto getSpaceSize(const SharedSpace<SpaceType>& s) {
    return getSpaceSize(*s);
  }

  template <QuadratureSpaceConcept SpaceType>
  auto getNumberOfCells(const SharedSpace<SpaceType>& s) {
    return getNumberOfCells(*s);
  }

  template <QuadratureSpaceConcept SpaceType>
  auto getNumberOfQuadraturePoints(const SharedSpace<SpaceType>& s,
                                   const cell_index<SpaceType>& e) {
    return getNumberOfQuadraturePoints(*s, e);
  }

  template <LinearQuadratureSpaceConcept SpaceType>
  auto getQuadraturePointOffset(const SharedSpace<SpaceType>& s,
                                const cell_index<SpaceType>& e,
                                const quadrature_point_index<SpaceType>& i) {
    return getQuadraturePointOffset(*s, e, i);
  }

}  // end of namespace mgis::function

#endif /* LIB_MFEM_MGIS_SHAREDSPACE_HXX */
