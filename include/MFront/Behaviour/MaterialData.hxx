/*!
 * \file   MaterialData.hxx
 * \brief
 * \author th202608
 * \date   19/06/2018
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
