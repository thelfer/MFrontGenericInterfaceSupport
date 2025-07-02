/*!
 * \file   MGIS/Behaviour/MaterialFunctionManager.ixx
 * \brief    
 * \author Thomas Helfer
 * \date   02/07/2025
 */

#ifndef LIB_MGIS_BEHAVIOUR_MATERIALFUNCTIONMANAGER_IXX
#define LIB_MGIS_BEHAVIOUR_MATERIALFUNCTIONMANAGER_IXX

namespace mgis::behaviour {

  template <
      mgis::function::LinearQuadratureSpaceConcept LinearQuadratureSpaceType>
  bool MaterialFunctionManager<LinearQuadratureSpaceType>::checkPreconditions(
      AbstractErrorHandler& eh,
      const std::shared_ptr<const LinearQuadratureSpaceType>& s,
      const Behaviour&) {
    if (s.get() == nullptr) {
      return eh.registerErrorMessage("invalid quadrature space");
    }
    return true;
  }  // end of checkPreconditions

  template <
      mgis::function::LinearQuadratureSpaceConcept LinearQuadratureSpaceType>
  MaterialFunctionManager<LinearQuadratureSpaceType>::MaterialFunctionManager(
      const std::shared_ptr<const LinearQuadratureSpaceType>& s,
      const Behaviour& behaviour)
      : MaterialFunctionManager(preconditions_check, s, behaviour) {}

  template <
      mgis::function::LinearQuadratureSpaceConcept LinearQuadratureSpaceType>
  template <bool doPreconditionsCheck>
  MaterialFunctionManager<LinearQuadratureSpaceType>::MaterialFunctionManager(
      const PreconditionsCheck<doPreconditionsCheck>& pcheck,
      const std::shared_ptr<const LinearQuadratureSpaceType>& s,
      const Behaviour& behaviour)
      : PreconditionsChecker<MaterialFunctionManager>(pcheck, s, behaviour),
        MaterialDataManager(static_cast<size_type>(getSpaceSize(s)), behaviour),
        qspace(s) {}  // end of MaterialFunctionManager

}  // end of namespace mgis::behaviour

#endif /* LIB_MGIS_BEHAVIOUR_MATERIALFUNCTIONMANAGER_IXX */
