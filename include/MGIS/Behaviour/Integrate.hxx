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
     * \brief structure in charge of handling temporary memory access.
     */
    struct MGIS_EXPORT IntegrateWorkSpace {
      /*!
       * \brief constructor
       * \param[in] b: behaviour
       */
      IntegrateWorkSpace(const Behaviour&);
      IntegrateWorkSpace(IntegrateWorkSpace&&);
      IntegrateWorkSpace(const IntegrateWorkSpace&);
      IntegrateWorkSpace& operator=(IntegrateWorkSpace&&);
      IntegrateWorkSpace& operator=(const IntegrateWorkSpace&);
      //! material properties at the beginning of the time step
      std::vector<real> mps0;
      //! material properties at the end of the time step
      std::vector<real> mps1;
      //! external state variables at the beginning of the time step
      std::vector<real> esvs0;
      //! external state variables at the end of the time step
      std::vector<real> esvs1;
    };  // end of struct IntegrateWorkSpace

    /*!
     * \brief return a thread-specific workspace associated with the given
     * behaviour.
     * \param[in] b: behaviour
     */
    MGIS_EXPORT IntegrateWorkSpace& getIntegrateWorkSpace(const Behaviour&);

    /*!
     * \brief integrate the behaviour. The returned value has the following
     * meaning:
     * - -1: integration failed
     * -  0: integration succeeded but results are unreliable
     * -  1: integration succeeded and results are reliable
     *
     * \param[in,out] d: behaviour data
     * \param[in,out] b: behaviour
     */
    int integrate(BehaviourDataView&, const Behaviour&);

    /*!
     * \brief integrate the behaviour for a range of integration points. The
     * returned value has the following meaning:
     * - -1: integration failed for at least one Gauss point
     * -  0: all integrations succeeded but results are unreliable for at least
     *       one Gauss point
     * -  1: integration succeeded and results are reliable.
     *
     * \param[in,out] m: material data manager
     * \param[in] dt: time step
     * \param[in] b: first index of the range
     * \param[in] e: last index of the range
     */
    MGIS_EXPORT int integrate(MaterialDataManager&,
                              const real,
                              const size_type,
                              const size_type);

    /*!
     * \brief integrate the behaviour for a range of integration points. The
     * returned value has the following meaning:
     * - -1: integration failed for at least one Gauss point
     * -  0: all integrations succeeded but results are unreliable for at least
     *       one Gauss point
     * -  1: integration succeeded and results are reliable.
     *
     * \param[in,out] p: thread pool
     * \param[in,out] m: material data manager
     * \param[in] dt: time step
     */
    MGIS_EXPORT int integrate(mgis::ThreadPool&,
                              MaterialDataManager&,
                              const real);

  }  // end of namespace behaviour

}  // end of namespace mgis

#include "MGIS/Behaviour/Integrate.ixx"

#endif /* LIB_MGIS_BEHAVIOUR_INTEGRATE_HXX */
