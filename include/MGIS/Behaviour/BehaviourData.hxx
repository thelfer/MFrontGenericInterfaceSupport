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

#include <iosfwd>
#include <vector>
#include "MGIS/Config.hxx"
#include "MGIS/Behaviour/State.hxx"
#include "MGIS/Behaviour/BehaviourDataView.hxx"

namespace mgis::behaviour {

  // forward declaration
  struct Behaviour;

  /*!
   * \brief structure in charge of containing the data required for a
   * behaviour integration.
   */
  struct MGIS_EXPORT BehaviourData {
    /*!
     * \brief constructor from a behaviour
     * \param[in] b: behaviour
     */
    BehaviourData(const Behaviour&);
    //! \brief move constructor
    BehaviourData(BehaviourData&&);
    //! \brief copy constructor
    BehaviourData(const BehaviourData&);
    //! \brief move assignement
    BehaviourData& operator=(BehaviourData&&);
    //! \brief copy assignement
    BehaviourData& operator=(const BehaviourData&);
    /*!
     * \brief a pointer to a buffer used to store error message
     *
     * By default, this is initialised to an internal buffer. This is
     * **not** thread-safe.
     *
     * The user may use its own buffer by setting this pointer appropriately
     * or set it to the null pointer. If not null, the pointer must
     * point to a buffer which is at least 512 characters wide (longer
     * error message are truncated).  The user must ensure
     * thread-safety (i.e. each thread shall have its own buffer).
     */
    char* error_message;
    //! \brief time increment
    mgis::real dt;
    /*!
     * \brief the stiffness matrix.
     *
     * On input, the first element of K (K[0]) must contain the type of type
     * of computation to be performed.
     *
     * Let Ke be equal to:
     *
     * - K[0] - 100 if K[0] is greater than 50
     * - K[0] otherwise.
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
    std::vector<real> K;
    /*!
     * \brief proposed time step increment increase factor
     *
     * The calling solver shall set a suitable value on input
     * depending on its policy befor each call to integrate.
     *
     * For instance, if the solver want to limit the increase to 20% at most, it
     * shall set it to 1.2. But setting it to 1, the solver won't allow the
     * behaviour to request an increase of the time step.
     */
    real rdt = 1;
    //! \brief speed of sound (only computed if requested)
    real speed_of_sound = 0;
    //! \brief state at the beginning of the time step
    State s0;
    //! \brief state at the end of the time step
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
  /*!
   * \brief print a detailled (verbose) description of the data associated
   * with an integration point using a markdown format
   * \param[in] os: ouptut stream
   * \param[in] b: behaviour
   * \param[in] d: behaviour data
   * \param[in] l: title level
   */
  MGIS_EXPORT void print_markdown(std::ostream&,
                                  const Behaviour&,
                                  const BehaviourData&,
                                  const mgis::size_type);

}  // end of namespace mgis::behaviour

#endif /* LIB_MGIS_BEHAVIOUR_BEHAVIOURDATA_HXX */
