/*!
 * \file   MGIS/Behaviour/MaterialFunctionManager.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   02/07/2025
 */

#ifndef LIB_MGIS_BEHAVIOUR_MATERIALFUNCTIONMANAGER_HXX
#define LIB_MGIS_BEHAVIOUR_MATERIALFUNCTIONMANAGER_HXX

#include <memory>
#include "MGIS/Config.hxx"
#include "MGIS/Contract.hxx"
#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"

namespace mgis::behaviour{

  /*!
   * \brief a simple wrapper around the MaterialDataManager class
   */
  template <
      mgis::function::LinearQuadratureSpaceConcept LinearQuadratureSpaceType>
  struct MaterialFunctionManager
      : public MaterialDataManager,
        private PreconditionsChecker<
            MaterialFunctionManager<LinearQuadratureSpaceType>> {
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] eh: error handler.
     * \param[in] s: quadrature space
     * \param[in] behaviour: behaviour.
     */
    [[nodiscard]] static bool checkPreconditions(
        AbstractErrorHandler&,
        const std::shared_ptr<const LinearQuadratureSpaceType>&,
        const Behaviour&);
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] behaviour: behaviour.
     */
    MaterialFunctionManager(
        const std::shared_ptr<const LinearQuadratureSpaceType>&,
        const Behaviour&);
    /*!
     * \brief constructor
     * \param[in] pcheck: do preconditions checks
     * \param[in] s: quadrature space.
     * \param[in] behaviour: behaviour.
     */
    template <bool doPreconditionsCheck>
    MaterialFunctionManager(
        const PreconditionsCheck<doPreconditionsCheck>&,
        const std::shared_ptr<const LinearQuadratureSpaceType>&,
        const Behaviour&);
    //! \brief return the quadrature space
    const LinearQuadratureSpaceType& getQuadratureSpace() const noexcept;
    //! \brief return the pointer to the quadrature space
    std::shared_ptr<const LinearQuadratureSpaceType> getQuadratureSpacePointer()
        const noexcept;

   private:
    //! \brief quadrature space
    std::shared_ptr<const LinearQuadratureSpaceType> qspace;
  };

} // end of namespace mgis::behaviour

#include "MGIS/Behaviour/MaterialFunctionManager.ixx"

#endif /* LIB_MGIS_BEHAVIOUR_MATERIALFUNCTIONMANAGER_HXX */
