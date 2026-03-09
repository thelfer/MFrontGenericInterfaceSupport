/*!
 * \file   MGIS/Behaviour/TimeStepStage.hxx
 * \brief
 * \author Thomas Helfer
 * \date   24/01/2026
 */

#ifndef LIB_MGIS_BEHAVIOUR_TIMESTEPSTAGE_HXX
#define LIB_MGIS_BEHAVIOUR_TIMESTEPSTAGE_HXX

namespace mgis::behaviour {

  /*!
   * \brief an enumeration describing either the beginning of the time step or
   * the end of the time step
   */
  enum struct TimeStepStage { BEGINNING_OF_TIME_STEP, END_OF_TIME_STEP };
  //! \brief variable associated with the beginning of the time step
  inline constexpr auto bts = TimeStepStage::BEGINNING_OF_TIME_STEP;
  //! \brief variable associated with the end of the time step
  inline constexpr auto ets = TimeStepStage::END_OF_TIME_STEP;

}  // end of namespace mgis::behaviour

#endif /* LIB_MGIS_BEHAVIOUR_TIMESTEPSTAGE_HXX */
