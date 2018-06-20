/*!
 * \file   MFront/Behaviour/Variable.hxx
 * \brief
 * \author th202608
 * \date   19/06/2018
 */

#ifndef LIB_MFRONT_BEHAVIOUR_VARIABLE_HXX
#define LIB_MFRONT_BEHAVIOUR_VARIABLE_HXX

#include <string>
#include <vector>
#include "MFront/Config.hxx"
#include "MFront/Behaviour/Hypothesis.hxx"

namespace mfront {

  namespace behaviour {

    /*!
     * \brief structure describing a variable
     */
    struct Variable {
      //! name of the variable
      std::string name;
      //! type of the variable
      enum { SCALAR, STENSOR, TENSOR } type;
    };  // end of struct Description

    /*!
     * \return the size of an array that may contain the values described by the
     * given array of variables
     * \param[in] vs: variables
     * \param[in] h: modelling hypothesis
     */
    MFRONT_EXPORT size_type getArraySize(const std::vector<Variable>&,
                                         const Hypothesis);

  }  // end of namespace behaviour

}  // end of namespace mfront

#endif /* LIB_MFRONT_BEHAVIOUR_VARIABLE_HXX */
