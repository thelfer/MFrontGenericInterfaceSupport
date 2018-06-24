/*!
 * \file   MGIS/Behaviour/Variable.hxx
 * \brief
 * \author Thomas Helfer
 * \date   19/06/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_BEHAVIOUR_VARIABLE_HXX
#define LIB_MGIS_BEHAVIOUR_VARIABLE_HXX

#include <string>
#include <vector>
#include "MGIS/Config.hxx"
#include "MGIS/Behaviour/Hypothesis.hxx"

namespace mgis {

  namespace behaviour {

    /*!
     * \brief structure describing a variable
     */
    struct Variable {
      //! name of the variable
      std::string name;
      //! type of the variable
      enum { SCALAR, VECTOR, STENSOR, TENSOR } type;
    };  // end of struct Description

    /*!
     * \return the size of an array that may contain the values described by the
     * given array of variables
     * \param[in] vs: variables
     * \param[in] h: modelling hypothesis
     */
    MGIS_EXPORT size_type getArraySize(const std::vector<Variable>&,
                                         const Hypothesis);

  }  // end of namespace behaviour

}  // end of namespace mgis

#endif /* LIB_MGIS_BEHAVIOUR_VARIABLE_HXX */
