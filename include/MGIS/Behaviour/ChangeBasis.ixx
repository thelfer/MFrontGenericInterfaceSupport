/*!
 * \file   include/MGIS/Behaviour/ChangeBasis.ixx
 * \brief
 * \author Thomas Helfer
 * \date   03/09/2018
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_BEHAVIOUR_CHANGEBASIS_IXX
#define LIB_MGIS_BEHAVIOUR_CHANGEBASIS_IXX

namespace mgis::behaviour {

  template <typename Rotation>
  void changeBasis(real* const o,
                   const real* const i,
                   const std::vector<Variable>& vs,
                   const Hypothesis h,
                   const Rotation& r) {
    auto offset = size_type{};
    for (const auto& v : vs) {
      if (v.type == Variable::SCALAR) {
        ++offset;
      } else if (v.type == Variable::VECTOR) {
        r.rotateVector(o + offset, i + offset);
        offset += getSpaceDimension(h);
      } else if (v.type == Variable::STENSOR) {
        r.rotateStensor(o + offset, i + offset);
        offset += getStensorSize(h);
      } else if (v.type == Variable::TENSOR) {
        r.rotateTensor(o + offset, i + offset);
        offset += getTensorSize(h);
      }
    }
  }  // end of changeBasis

}  // end of namespace mgis::behaviour

#endif /* LIB_MGIS_BEHAVIOUR_CHANGEBASIS_IXX */
