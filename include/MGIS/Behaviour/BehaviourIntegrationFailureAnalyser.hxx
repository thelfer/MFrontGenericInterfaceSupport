/*!
 * \file   MGIS/Behaviour/BehaviourIntegrationFailureAnalyser.hxx
 * \brief  The header define some useful features for debugging purposes
 * \author Thomas Helfer
 * \date   21/08/2025
 */

#ifndef LIB_MGIS_BEHAVIOUR_BEHAVIOURINTEGRATIONFAILUREANALYSER_HXX
#define LIB_MGIS_BEHAVIOUR_BEHAVIOURINTEGRATIONFAILUREANALYSER_HXX

#include <string>
#include <functional>
#include <string_view>
#include "MGIS/Config.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"
#include "MGIS/Behaviour/BehaviourDataView.hxx"

namespace mgis::behaviour {

  // forward declarations
  struct Behaviour;

}  // end of namespace mgis::behaviour

namespace mgis::behaviour::debug {

  /*!
   * \brief object in charge of analysing and integration failure. This object
   * is used by the `integrate_debug` function
   *
   * \note most analysers generate a file. See the
   * `getBehaviourIntegrationFailureAnalysisFileName` below for a way to
   * generate a unique file name. By default, this file name is unique for the
   * given process. In an MPI context, one shall provide a file name generator
   * which takes the MPI process id into account to avoid conflicts, see the
   * function `setBehaviourIntegrationFailureAnalysisFileNameGenerator` for
   * details.
   */
  struct MGIS_EXPORT BehaviourIntegrationFailureAnalyser {
    /*!
     * \brief handle an integration failure
     * \param[in] b: behaviour being integrated
     * \param[in] d: behaviour data leading to the integration failure
     */
    virtual void analyse(const Behaviour&,
                         const BehaviourData&) const noexcept = 0;
    /*!
     * \brief handle an integration failure
     * \param[in] b: behaviour being integrated
     * \param[in] d: behaviour data view leading to the integration failure
     */
    virtual void analyse(const Behaviour&,
                         const BehaviourDataView&) const noexcept = 0;
    /*!
     * \brief use a copy of the inputs
     *
     * If true, the behaviour data view  passed to the `integrate_debug`
     * function is copied before integrating the behaviour. Copying has
     * a significant overhead and is normally not required for correctly written
     * behaviours but data corruption may occur, so copying is safer.
     */
    virtual bool shallCopyBehaviourDataBeforeIntegration() const noexcept = 0;
    //! \brief destructor
    virtual ~BehaviourIntegrationFailureAnalyser() noexcept;
  };

  /*!
   * \return a file name used to output
   * debugging information in case of integration
   * failure
   *
   * \param[in] n: behaviour name
   * \param[in] ext: file name
   */
  MGIS_EXPORT std::string getBehaviourIntegrationFailureAnalysisFileName(
      std::string_view, std::string_view);

  /*!
   * \brief a function hook used to generate a file name
   * for debugging behaviour integration failures
   * \param[in] g: generator
   *
   * \note the generator takes the behaviour name as first argument and a unique
   * identifier (for the given process) as second argument
   */
  MGIS_EXPORT void setBehaviourIntegrationFailureAnalysisFileNameGenerator(
      std::function<std::string(std::string_view,  // behaviour name
                                std::string_view,  // unique identifier
                                std::string_view)  // file extension
                    >);

#ifndef LIB_MGIS_BEHAVIOUR_INTEGRATE_HXX
  /*!
   * \return the default debugging options
   */
  MGIS_EXPORT const debug::BehaviourIntegrationFailureAnalyser&
  getDefaultBehaviourIntegrationFailureAnalyser();
#endif

}  // end of namespace mgis::behaviour::debug

#endif /* LIB_MGIS_BEHAVIOUR_BEHAVIOURINTEGRATIONFAILUREANALYSER_HXX */
