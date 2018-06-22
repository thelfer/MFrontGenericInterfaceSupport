/*!
 * \file   MaterialData.hxx
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

#ifndef LIB_MFRONT_BEHAVIOUR_DATA_HXX
#define LIB_MFRONT_BEHAVIOUR_DATA_HXX

#include "MFront/Behaviour/Information.hxx"

namespace mfront {

  /*!
   * \brief structure containing the material data of a set of integration points
   * \tparam real: numerical type used
   */
  struct MaterialData {
    //! \brief material properties
    std::vector<real> mps;
    //! \brief internal state variables
    std::vector<real> ivs;
    //! \brief external state variables
    std::vector<real> evs;
  }; // end of struct MaterialData

}  // end of mfront

#endif /* LIB_MFRONT_BEHAVIOUR_DATA_HXX */
