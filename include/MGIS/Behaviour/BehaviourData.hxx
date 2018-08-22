/*!
 * \file   BehaviourData.hxx
 * \brief
 * \author Thomas Helfer
 * \date   02/08/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_BEHAVIOUR_BEHAVIOURDATA_HXX
#define LIB_MGIS_BEHAVIOUR_BEHAVIOURDATA_HXX

#include <vector>
#include "MGIS/Config.hxx"
#include "MGIS/Behaviour/State.hxx"
#include "MGIS/Behaviour/BehaviourDataView.hxx"

namespace mgis {

  namespace behaviour {

    // forward declaration
    struct Behaviour;

    /*!
     * \brief structure in charge of containing the data required for a
     * behaviour integration:
     */
    struct MGIS_EXPORT BehaviourData {
      /*!
       * \brief constructor from a behaviour
       * \param[in] b: behaviour
       */
      BehaviourData(const Behaviour&);
      //! move constructor
      BehaviourData(BehaviourData&&);
      //! copy constructor
      BehaviourData(const BehaviourData&);
      //! move assignement
      BehaviourData& operator=(BehaviourData&&);
      //! copy assignement
      BehaviourData& operator=(const BehaviourData&);
      //! time increment
      real dt;
      /*!
       * \brief the stiffness matrix.
       * On input, the first element of K must contain the type of type
       * of stiffness matrix expected. If this value is negative, only
       * the prediction operator is computed. This value has the
       * following meaning:
       * - if K[0] is lower than -2.5, the tangent operator must be
       *   computed.
       * - if K[0] is in [-2.5:-1.5]: the secant operator must be
       *   computed.
       * - if K[0] is in [-1.5:-0.5]: the elastic operator must be
       *   computed.
       * - if K[0] is in [-0.5:0.5]: the behaviour integration is
       *   performed, but no stiffness matrix.
       * - if K[0] is in [0.5:1.5]: the elastic operator must be
       *   computed.
       * - if K[0] is in [1.5:2.5]: the secant operator must be
       *   computed.
       * - if K[0] is in [2.5:3.5]: the secant operator must be
       *   computed.
       * - if K[0] is greater than 3.5, the consistent tangent operator
       *   must be computed.
       */
      std::vector<real> K;
      /*!
       * \brief proposed time step increment increase factor
       */
      real rdt;
      /*!
       * \brief state at the beginning of the time step
       */
      State s0;
      /*!
       * \brief state at the end of the time step
       */
      State s1;
    };  // end of struct BehaviourData

    /*!
     * \brief update the behaviour data by:
     * - setting s1 equal to s0
     * - filling the stiffness matrix with 0
     */
    MGIS_EXPORT void update(BehaviourData&);
    /*!
     * \brief revert the behaviour data by:
     * - setting s1 equal to s0
     * - filling the stiffness matrix with 0
     */
    MGIS_EXPORT void revert(BehaviourData&);

    /*!
     * \brief make a view from a behaviour data
     * \param[in] d: data
     * \return the view
     * \note the view has at most the same life time as the data.
     */
    MGIS_EXPORT BehaviourDataView make_view(BehaviourData&);

  }  // end of namespace behaviour

}  // end of namespace mgis

#endif /* LIB_MGIS_BEHAVIOUR_BEHAVIOURDATA_HXX */
