/*!
 * \file   MGIS/Behaviour/MaterialFunctionManager.hxx
 * \brief
 * \author Thomas Helfer
 * \date   02/07/2025
 */

#ifndef LIB_MGIS_BEHAVIOUR_MATERIALFUNCTIONMANAGER_HXX
#define LIB_MGIS_BEHAVIOUR_MATERIALFUNCTIONMANAGER_HXX

#include <memory>
#include <optional>
#include <string_view>
#include "MGIS/Config.hxx"
#include "MGIS/Contract.hxx"
#include "MGIS/Function/SpaceConcept.hxx"
#include "MGIS/Function/SharedSpace.hxx"
#include "MGIS/Function/Function.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"

namespace mgis::behaviour {

  /*!
   * \brief a simple wrapper around the MaterialDataManager class
   */
  template <mgis::function::LinearElementSpaceConcept SpaceType>
  struct MaterialFunctionManager
      : private PreconditionsChecker<MaterialFunctionManager<SpaceType>>,
        public MaterialDataManager {
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] eh: error handler.
     * \param[in] s: quadrature space
     * \param[in] behaviour: behaviour.
     */
    [[nodiscard]] static bool checkPreconditions(
        AbstractErrorHandler&,
        const std::shared_ptr<const SpaceType>&,
        const Behaviour&);
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] behaviour: behaviour.
     */
    MaterialFunctionManager(const std::shared_ptr<const SpaceType>&,
                            const Behaviour&);
    /*!
     * \brief constructor
     * \param[in] pcheck: do preconditions checks
     * \param[in] s: quadrature space.
     * \param[in] behaviour: behaviour.
     */
    template <bool doPreconditionsCheck>
    MaterialFunctionManager(const PreconditionsCheck<doPreconditionsCheck>&,
                            const std::shared_ptr<const SpaceType>&,
                            const Behaviour&);
    /*!
     * \brief check that the preconditions to build the view are met
     * \param[in] eh: error handler.
     * \param[in] s: quadrature space
     * \param[in] behaviour: behaviour.
     */
    [[nodiscard]] static bool checkPreconditions(
        AbstractErrorHandler&,
        const std::shared_ptr<const SpaceType>&,
        const std::shared_ptr<const Behaviour>&);
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] behaviour: behaviour.
     */
    MaterialFunctionManager(const std::shared_ptr<const SpaceType>&,
                            const std::shared_ptr<const Behaviour>&);
    /*!
     * \brief constructor
     * \param[in] pcheck: do preconditions checks
     * \param[in] s: quadrature space.
     * \param[in] behaviour: behaviour.
     */
    template <bool doPreconditionsCheck>
    MaterialFunctionManager(const PreconditionsCheck<doPreconditionsCheck>&,
                            const std::shared_ptr<const SpaceType>&,
                            const std::shared_ptr<const Behaviour>&);
    //! \brief return the quadrature space
    const SpaceType& getSpace() const noexcept;
    //! \brief return the pointer to the quadrature space
    std::shared_ptr<const SpaceType> getSpacePointer() const noexcept;
    //! \brief return the pointer to the quadrature space
    mgis::function::SharedSpace<SpaceType> getSharedSpace() const noexcept;

   private:
    //! \brief quadrature space
    std::shared_ptr<const SpaceType> qspace;
    //! \brief behaviour, may be null if the behaviour is not managed
    std::shared_ptr<const Behaviour> behaviour_ptr;
  };

  enum TimeStepStage {
    BEGINNING_OF_TIME_STEP,
    END_OF_TIME_STEP,
  };

  inline constexpr auto bts = TimeStepStage::BEGINNING_OF_TIME_STEP;
  inline constexpr auto ets = TimeStepStage::END_OF_TIME_STEP;

  // this variable is a workaround Visual Studio 2022 limitation
  // on defining default value for non-type template argument
  template <size_type N>
  inline constexpr auto fixed_size_dynamic_stride_data_layout_description =
      ::mgis::function::FunctionDataLayoutDescription{
          .data_size = N, .data_stride = dynamic_extent};

  /*!
   * \brief return the gradient associated with the given name
   * \param[in] eh: error handler
   * \param[in] m: material function manager
   * \param[in] n: name of the gradient
   * \param[in] ts: time step stage
   */
  template <mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<
      mgis::function::FunctionView<mgis::function::SharedSpace<SpaceType>>>
  getGradient(AbstractErrorHandler&,
              MaterialFunctionManager<SpaceType>&,
              std::string_view,
              const TimeStepStage = ets);
  /*!
   * \brief return the gradient associated with the given name
   * \param[in] eh: error handler
   * \param[in] m: material function manager
   * \param[in] n: name of the gradient
   * \param[in] ts: time step stage
   */
  template <mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<
      mgis::function::
          FunctionView<mgis::function::SharedSpace<SpaceType>, {}, false>>
  getGradient(AbstractErrorHandler&,
              const MaterialFunctionManager<SpaceType>&,
              std::string_view,
              const TimeStepStage = ets);
  /*!
   * \brief return the gradient associated with the given name
   * \param[in] eh: error handler
   * \param[in] m: material function manager
   * \param[in] n: name of the gradient
   * \param[in] ts: time step stage
   */
  template <size_type N, mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<mgis::function::FunctionView<
      mgis::function::SharedSpace<SpaceType>,
      fixed_size_dynamic_stride_data_layout_description<N>>>
  getGradient(AbstractErrorHandler&,
              MaterialFunctionManager<SpaceType>&,
              std::string_view,
              const TimeStepStage = ets);
  /*!
   * \brief return the gradient associated with the given name
   * \param[in] eh: error handler
   * \param[in] m: material function manager
   * \param[in] n: name of the gradient
   * \param[in] ts: time step stage
   */
  template <size_type N, mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<mgis::function::FunctionView<
      mgis::function::SharedSpace<SpaceType>,
      fixed_size_dynamic_stride_data_layout_description<N>,
      false>>
  getGradient(AbstractErrorHandler&,
              const MaterialFunctionManager<SpaceType>&,
              std::string_view,
              const TimeStepStage = ets);

  /*!
   * \brief return the thermodynamic force associated with the given name
   * \param[in] eh: error handler
   * \param[in] m: material function manager
   * \param[in] n: name of the thermodynamic force
   * \param[in] ts: time step stage
   */
  template <mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<
      mgis::function::FunctionView<mgis::function::SharedSpace<SpaceType>>>
  getThermodynamicForce(AbstractErrorHandler&,
                        MaterialFunctionManager<SpaceType>&,
                        std::string_view,
                        const TimeStepStage = ets);
  /*!
   * \brief return the thermodynamic force associated with the given name
   * \param[in] eh: error handler
   * \param[in] m: material function manager
   * \param[in] n: name of the thermodynamic force
   * \param[in] ts: time step stage
   */
  template <mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<
      mgis::function::
          FunctionView<mgis::function::SharedSpace<SpaceType>, {}, false>>
  getThermodynamicForce(AbstractErrorHandler&,
                        const MaterialFunctionManager<SpaceType>&,
                        std::string_view,
                        const TimeStepStage = ets);
  /*!
   * \brief return the thermodynamic force associated with the given name
   * \param[in] eh: error handler
   * \param[in] m: material function manager
   * \param[in] n: name of the thermodynamic force
   * \param[in] ts: time step stage
   */
  template <size_type N, mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<mgis::function::FunctionView<
      mgis::function::SharedSpace<SpaceType>,
      fixed_size_dynamic_stride_data_layout_description<N>>>
  getThermodynamicForce(AbstractErrorHandler&,
                        MaterialFunctionManager<SpaceType>&,
                        std::string_view,
                        const TimeStepStage = ets);
  /*!
   * \brief return the thermodynamic force associated with the given name
   * \param[in] eh: error handler
   * \param[in] m: material function manager
   * \param[in] n: name of the thermodynamic force
   * \param[in] ts: time step stage
   */
  template <size_type N, mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<mgis::function::FunctionView<
      mgis::function::SharedSpace<SpaceType>,
      fixed_size_dynamic_stride_data_layout_description<N>,
      false>>
  getThermodynamicForce(AbstractErrorHandler&,
                        const MaterialFunctionManager<SpaceType>&,
                        std::string_view,
                        const TimeStepStage = ets);

  /*!
   * \brief return the internal state variable associated with the given name
   * \param[in] eh: error handler
   * \param[in] m: material function manager
   * \param[in] n: name of the internal state variable
   * \param[in] ts: time step stage
   */
  template <mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<
      mgis::function::FunctionView<mgis::function::SharedSpace<SpaceType>>>
  getInternalStateVariable(AbstractErrorHandler&,
                           MaterialFunctionManager<SpaceType>&,
                           std::string_view,
                           const TimeStepStage = ets);
  /*!
   * \brief return the internal state variable associated with the given name
   * \param[in] eh: error handler
   * \param[in] m: material function manager
   * \param[in] n: name of the internal state variable
   * \param[in] ts: time step stage
   */
  template <mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<
      mgis::function::
          FunctionView<mgis::function::SharedSpace<SpaceType>, {}, false>>
  getInternalStateVariable(AbstractErrorHandler&,
                           const MaterialFunctionManager<SpaceType>&,
                           std::string_view,
                           const TimeStepStage = ets);
  /*!
   * \brief return the internal state variable associated with the given name
   * \param[in] eh: error handler
   * \param[in] m: material function manager
   * \param[in] n: name of the internal state variable
   * \param[in] ts: time step stage
   */
  template <size_type N, mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<mgis::function::FunctionView<
      mgis::function::SharedSpace<SpaceType>,
      fixed_size_dynamic_stride_data_layout_description<N>>>
  getInternalStateVariable(AbstractErrorHandler&,
                           MaterialFunctionManager<SpaceType>&,
                           std::string_view,
                           const TimeStepStage = ets);
  /*!
   * \brief return the internal state variable associated with the given name
   * \param[in] eh: error handler
   * \param[in] m: material function manager
   * \param[in] n: name of the internal state variable
   * \param[in] ts: time step stage
   */
  template <size_type N, mgis::function::LinearElementSpaceConcept SpaceType>
  std::optional<mgis::function::FunctionView<
      mgis::function::SharedSpace<SpaceType>,
      fixed_size_dynamic_stride_data_layout_description<N>,
      false>>
  getInternalStateVariable(AbstractErrorHandler&,
                           const MaterialFunctionManager<SpaceType>&,
                           std::string_view,
                           const TimeStepStage = ets);

}  // end of namespace mgis::behaviour

#include "MGIS/Behaviour/MaterialFunctionManager.ixx"

#endif /* LIB_MGIS_BEHAVIOUR_MATERIALFUNCTIONMANAGER_HXX */
