/*!
 * \file   Hypothesis.hxx
 * \brief
 * \author th202608
 * \date   19/06/2018
 */

#ifndef LIB_MFRONT_BEHAVIOUR_HYPOTHESIS_HXX
#define LIB_MFRONT_BEHAVIOUR_HYPOTHESIS_HXX

#include <string>
#include "MFront/Config.hxx"

namespace mfront {

  namespace behaviour {

    //! \brief the list of supported modelling hypotheses
    enum struct Hypothesis {
      AXISYMMETRICALGENERALISEDPLANESTRAIN,
      AXISYMMETRICALGENERALISEDPLANESTRESS,
      AXISYMMETRICAL,
      PLANESTRESS,
      PLANESTRAIN,
      GENERALISEDPLANESTRAIN,
      TRIDIMENSIONAL
    };  // end of enum Hypothesis

    /*!
     * \return the string associated to the given hypothesis
     * \param[in] h: modelling hypothesis
     */
    MFRONT_EXPORT const char* toString(const Hypothesis);

  }  // end of namespace behaviour

}  // end of namespace mfront

#endif /* LIB_MFRONT_BEHAVIOUR_HYPOTHESIS_HXX */
