/*!
 * \file   Integrate.hxx
 * \brief
 * \author Thomas Helfer
 * \date   01/08/2018
 * \copyright Copyright (C) 2006-2018 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence
 * or the CECILL-A licence. A copy of thoses licences are delivered
 * with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MGIS_BEHAVIOUR_INTEGRATE_HXX
#define LIB_MGIS_BEHAVIOUR_INTEGRATE_HXX

#include <limits>
#include <thread>
#include <vector>
#include "MGIS/Config.hxx"
#include "MGIS/Behaviour/BehaviourDataView.hxx"

namespace mgis {

  // forward declaration
  struct ThreadPool;

  namespace behaviour {

    // forward declaration
    struct Behaviour;
    // forward declaration
    struct MaterialDataManager;

    /*!
     * \brief type of integration to be performed
     */
    enum struct IntegrationType {
      PREDICTION_TANGENT_OPERATOR = -3,
      PREDICTION_SECANT_OPERATOR = -2,
      PREDICTION_ELASTIC_OPERATOR = -1,
      INTEGRATION_NO_TANGENT_OPERATOR = 0,
      INTEGRATION_ELASTIC_OPERATOR = 1,
      INTEGRATION_SECANT_OPERATOR = 2,
      INTEGRATION_TANGENT_OPERATOR = 3,
      INTEGRATION_CONSISTENT_TANGENT_OPERATOR = 4
    };  // end of enum IntegrationType

    /*!
   * \brief structure defining various option
   */
  struct BehaviourIntegrationOptions {
    //! \brief type of integration to be performed
    IntegrationType integration_type =
        IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR;
    //! \brief if true, the speed of sound shall be computed
    bool compute_speed_of_sound = false;
  };  // end of BehaviourIntegrationOptions

  /*!
   * \brief structure in charge of reporting the result of a behaviour
   * integration.
   */
  struct MGIS_EXPORT BehaviourIntegrationResult {
    //! \brief default constructor
    BehaviourIntegrationResult();
    //! \brief move constructor
    BehaviourIntegrationResult(BehaviourIntegrationResult&&);
    //! \brief copye constructor
    BehaviourIntegrationResult(const BehaviourIntegrationResult&);
    //! \brief move assignement
    BehaviourIntegrationResult& operator=(BehaviourIntegrationResult&&);
    //! \brief copy assignement
    BehaviourIntegrationResult& operator=(const BehaviourIntegrationResult&);
    //! \brief destructor
    ~BehaviourIntegrationResult();
    /*! \brief exit status
     *
     * The returned value has the following meaning:
     * - -1: integration failed for at least one Gauss point
     * -  0: all integrations succeeded but results are unreliable for at least
     *       one Gauss point
     * -  1: integration succeeded and results are reliable.
     */
    int exit_status = 1;
    //! \brief proposed time step increase factor
    mgis::real time_step_increase_factor = 10;
    /*!
     * \brief number of the integration point that failed or number of
     * integration that reported unreliable results.
     */
    mgis::size_type n = std::numeric_limits<mgis::size_type>::max();
    //! \brief error message, if any
    std::string error_message;
  };  // end of struct BehaviourIntegrationResult

  /*!
   * \brief structure in charge of reporting the result of a behaviour
   * integration.
   */
  struct MGIS_EXPORT MultiThreadedBehaviourIntegrationResult {
    //! \brief default constructor
    MultiThreadedBehaviourIntegrationResult();
    //! \brief move constructor
    MultiThreadedBehaviourIntegrationResult(
        MultiThreadedBehaviourIntegrationResult&&);
    //! \brief copye constructor
    MultiThreadedBehaviourIntegrationResult(
        const MultiThreadedBehaviourIntegrationResult&);
    //! \brief move assignement
    MultiThreadedBehaviourIntegrationResult& operator=(
        MultiThreadedBehaviourIntegrationResult&&);
    //! \brief copy assignement
    MultiThreadedBehaviourIntegrationResult& operator=(
        const MultiThreadedBehaviourIntegrationResult&);
    //! \brief destructor
    ~MultiThreadedBehaviourIntegrationResult();
    /*! \brief exit status
     *
     * The returned value has the following meaning:
     * - -1: integration failed for at least one Gauss point
     * -  0: all integrations succeeded but results are unreliable for at least
     *       one Gauss point
     * -  1: integration succeeded and results are reliable.
     */
    int exit_status = 1;
    //! \brief integration results per threads
    std::vector<BehaviourIntegrationResult> results;
  };  // end of struct MultiThreadedBehaviourIntegrationResult

  /*!
   * \brief integrate the behaviour. The returned value has the following
   * meaning:
   * - -1: integration failed
   * -  0: integration succeeded but results are unreliable
   * -  1: integration succeeded and results are reliable
   *
   * \param[in,out] d: behaviour data
   * \param[in,out] b: behaviour
   *

   * \note: the type of integration to be performed, must be
   * explicitely set in d.K[0], as follows (see the `IntegrationType` enum).
   *
   * Let Ke be equal to:
   *
   * - d.K[0] - 100 if d.K[0] is greater than 50
   * - d.K[0] otherwise.
   *
   * If Ke is negative, only the prediction operator is computed and
   * no behaviour integration is performed.
   *
   * Ke has the following meaning:
   *
   * - if Ke is lower than -2.5, the tangent operator must be
   *   computed.
   * - if Ke is in [-2.5:-1.5]: the secant operator must be
   *   computed.
   * - if Ke is in [-1.5:-0.5]: the elastic operator must be
   *   computed.
   * - if Ke is in [-0.5:0.5]: the behaviour integration is
   *   performed, but no stiffness matrix.
   * - if Ke is in [0.5:1.5]: the elastic operator must be
   *   computed.
   * - if Ke is in [1.5:2.5]: the secant operator must be
   *   computed.
   * - if Ke is in [2.5:3.5]: the secant operator must be
   *   computed.
   * - if Ke is greater than 3.5, the consistent tangent operator
   *   must be computed.
   */
    int integrate(BehaviourDataView&, const Behaviour&);

    /*!
     * \brief integrate the behaviour for a range of integration points.
   * \return the result of the behaviour integration.
   * \param[in,out] m: material data manager
   * \param[in] c: description of the operation to be performed
   * \param[in] dt: time step
   * \param[in] b: first index of the range
   * \param[in] e: last index of the range
   *
   * \note if required, the memory associated with the tangent operator blocks
   * is automatically allocated.
   */
  MGIS_EXPORT BehaviourIntegrationResult
  integrate(MaterialDataManager&,
            const BehaviourIntegrationOptions&,
            const real,
            const size_type,
            const size_type);

  /*!
   * \brief integrate the behaviour over all integration points using a thread
   * pool to parallelize the integration.
   * \return the result of the behaviour integration.
   * \param[in,out] p: thread pool
   * \param[in,out] m: material data manager
   * \param[in] c: description of the operation to be performed
   * \param[in] dt: time step
   *
   * \note if required, the memory associated with the tangent operator blocks
   * is automatically allocated.
   */
  MGIS_EXPORT MultiThreadedBehaviourIntegrationResult
  integrate(mgis::ThreadPool&,
            MaterialDataManager&,
            const BehaviourIntegrationOptions&,
            const real);

  /*!
   * \brief integrate the behaviour for a range of integration points.
     * \return an exit status. The returned value has the following meaning:
     * - -1: integration failed for at least one Gauss point
     * -  0: all integrations succeeded but results are unreliable for at least
     *       one Gauss point
     * -  1: integration succeeded and results are reliable.
     *
     * \param[in,out] m: material data manager
     * \param[in] dt: time step
     * \param[in] b: first index of the range
     * \param[in] e: last index of the range
     *
     * \note if required, the memory associated with the tangent operator blocks
     * is automatically allocated.
     */
    MGIS_EXPORT int integrate(MaterialDataManager&,
                              const IntegrationType,
                              const real,
                              const size_type,
                              const size_type);

    /*!
     * \brief integrate the behaviour for a range of integration points.
     * \return an exit status. The returned value has the following meaning:
     * - -1: integration failed for at least one integration point
     * -  0: all integrations succeeded but results are unreliable for at least
     *       one Gauss point
     * -  1: integration succeeded and results are reliable.
     *
     * \param[in,out] p: thread pool
     * \param[in,out] m: material data manager
     * \param[in] dt: time step
     *
     * \note if required, the memory associated with the tangent operator blocks
     * is automatically allocated.
     */
    MGIS_EXPORT int integrate(mgis::ThreadPool&,
                              MaterialDataManager&,
                              const IntegrationType it,
                              const real);

  }  // end of namespace behaviour

}  // end of namespace mgis

#include "MGIS/Behaviour/Integrate.ixx"

#endif /* LIB_MGIS_BEHAVIOUR_INTEGRATE_HXX */
