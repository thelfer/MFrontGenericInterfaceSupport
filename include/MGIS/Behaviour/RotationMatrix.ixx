/*!
 * \file   include/MGIS/Behaviour/RotationMatrix.ixx
 * \brief
 * \author Thomas Helfer
 * \date   03/03/2021
 */

#ifndef LIB_MGIS_BEHAVIOUR_ROTATIONMATRIX_IXX
#define LIB_MGIS_BEHAVIOUR_ROTATIONMATRIX_IXX

namespace mgis::behaviour {

  inline std::array<mgis::real, 9u> buildRotationMatrix(
      const std::span<const mgis::real, 2u>& a) {
    return {a[0], -a[1], 0,  //
            a[1], a[0],  0,  //
            0,    0,     1};
  }  // end of buildRotationMatrix

  inline std::array<mgis::real, 9u> buildRotationMatrix(
      const std::span<const mgis::real, 3u>& a1,
      const std::span<const mgis::real, 3u>& a2) {
    return {a1[0], a2[0], a1[1] * a2[2] - a1[2] * a2[1],  //
            a1[1], a2[1], a1[2] * a2[0] - a1[0] * a2[2],  //
            a1[2], a2[2], a1[0] * a2[1] - a1[1] * a2[0]};
  }  // end of buildRotationMatrix

}  // end of namespace mgis::behaviour

#endif /* LIB_MGIS_BEHAVIOUR_ROTATIONMATRIX_IXX */
