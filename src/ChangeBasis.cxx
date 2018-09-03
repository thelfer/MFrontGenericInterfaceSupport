/*!
 * \file   ChangeBasis.cxx
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

#include "MGIS/Cste.hxx"
#include "MGIS/Behaviour/ChangeBasis.hxx"

    namespace mgis {

  namespace behaviour {

    Rotation2D::Rotation2D(const real* const m)
        : m00(m[0]),
          m01(m[1]),
          m10(m[3]),
          m11(m[4]) {}  // end of Rotation2D::Rotation2D
          
    void Rotation2D::rotateVector(real* const o, const real* const i) const {
      const real nv[2] = {
          m00 * i[0] + m01 * i[1], m10 * i[0] + m11 * i[1],
      };
      o[0] = nv[0];
      o[1] = nv[1];
    }  // end of Rotation2D::rotateVector

    void Rotation2D::rotateStensor(real* const o, const real* const i) const {
      constexpr const auto cste = mgis::Cste::sqrt2;
      const real nv[3] = {
          cste * m00 * m10 * i[3] + m10 * m10 * i[1] + m00 * m00 * i[0],
          cste * m01 * m11 * i[3] + m11 * m11 * i[1] + m01 * m01 * i[0],
          (m00 * m11 + m01 * m10) * i[3] +
              cste * (m10 * m11 * i[1] + m00 * m01 * i[0])};
      o[0] = nv[0];
      o[1] = nv[1];
      o[2] = i[2];
      o[3] = nv[2];
    }  // end of Rotation2D::rotateStensor

    void Rotation2D::rotateTensor(real* const o, const real* const i) const {
      const real nv[4u] = {m00 * m10 * i[4] + m00 * m10 * i[3] +
                               m10 * m10 * i[1] + m00 * m00 * i[0],
                           m01 * m11 * i[4] + m01 * m11 * i[3] +
                               m11 * m11 * i[1] + m01 * m01 * i[0],
                           m01 * m10 * i[4] + m00 * m11 * i[3] +
                               m10 * m11 * i[1] + m00 * m01 * i[0],
                           m00 * m11 * i[4] + m01 * m10 * i[3] +
                               m10 * m11 * i[1] + m00 * m01 * i[0]};
      o[0] = nv[0];
      o[1] = nv[1];
      o[2] = i[2];
      o[3] = nv[2];
      o[4] = nv[3];
    }  // end of Rotation2D::rotateTensor

    Rotation3D::Rotation3D(const real* const m)
        : m00(m[0]),
          m01(m[1]),
          m02(m[2]),
          m10(m[3]),
          m11(m[4]),
          m12(m[5]),
          m20(m[6]),
          m21(m[7]),
          m22(m[8]) {}  // end of Rotation3D::Rotation3D

    void Rotation3D::rotateVector(real* const o, const real* const i) const {
      const real nv[3] = {
          m00 * i[0] + m01 * i[1] + m02 * i[2],
          m10 * i[0] + m11 * i[1] + m12 * i[2],
          m20 * i[0] + m21 * i[1] + m22 * i[2],
      };
      o[0] = nv[0];
      o[1] = nv[1];
      o[2] = nv[2];
    }  // end of Rotation3D::rotateVector

    void Rotation3D::rotateStensor(real* const o, const real* const i) const {
      constexpr const auto cste = mgis::Cste::sqrt2;
      const real nv[6] = {
          cste * (m10 * m20 * i[5] + m00 * m20 * i[4] + m00 * m10 * i[3]) +
              m20 * m20 * i[2] + m10 * m10 * i[1] + m00 * m00 * i[0],
          cste * (m11 * m21 * i[5] + m01 * m21 * i[4] + m01 * m11 * i[3]) +
              m21 * m21 * i[2] + m11 * m11 * i[1] + m01 * m01 * i[0],
          cste * (m12 * m22 * i[5] + m02 * m22 * i[4] + m02 * m12 * i[3]) +
              m22 * m22 * i[2] + m12 * m12 * i[1] + m02 * m02 * i[0],
          (m10 * m21 + m11 * m20) * i[5] + (m00 * m21 + m01 * m20) * i[4] +
              (m00 * m11 + m01 * m10) * i[3] +
              cste * (m20 * m21 * i[2] + m10 * m11 * i[1] + m00 * m01 * i[0]),
          (m10 * m22 + m12 * m20) * i[5] + (m00 * m22 + m02 * m20) * i[4] +
              (m00 * m12 + m02 * m10) * i[3] +
              cste * (m20 * m22 * i[2] + m10 * m12 * i[1] + m00 * m02 * i[0]),
          (m11 * m22 + m12 * m21) * i[5] + (m01 * m22 + m02 * m21) * i[4] +
              (m01 * m12 + m02 * m11) * i[3] +
              cste * (m21 * m22 * i[2] + m11 * m12 * i[1] + m01 * m02 * i[0])};
      o[0] = nv[0];
      o[1] = nv[1];
      o[2] = nv[2];
      o[3] = nv[3];
      o[4] = nv[4];
      o[5] = nv[5];
    }  // end of Rotation3D::rotateStensor

    void Rotation3D::rotateTensor(real* const o, const real* const i) const {
      const real nv[9] = {
          m10 * m20 * i[8] + m10 * m20 * i[7] + m00 * m20 * i[6] +
              m00 * m20 * i[5] + m00 * m10 * i[4] + m00 * m10 * i[3] +
              m20 * m20 * i[2] + m10 * m10 * i[1] + m00 * m00 * i[0],
          m11 * m21 * i[8] + m11 * m21 * i[7] + m01 * m21 * i[6] +
              m01 * m21 * i[5] + m01 * m11 * i[4] + m01 * m11 * i[3] +
              m21 * m21 * i[2] + m11 * m11 * i[1] + m01 * m01 * i[0],
          m12 * m22 * i[8] + m12 * m22 * i[7] + m02 * m22 * i[6] +
              m02 * m22 * i[5] + m02 * m12 * i[4] + m02 * m12 * i[3] +
              m22 * m22 * i[2] + m12 * m12 * i[1] + m02 * m02 * i[0],
          m11 * m20 * i[8] + m10 * m21 * i[7] + m01 * m20 * i[6] +
              m00 * m21 * i[5] + m01 * m10 * i[4] + m00 * m11 * i[3] +
              m20 * m21 * i[2] + m10 * m11 * i[1] + m00 * m01 * i[0],
          m10 * m21 * i[8] + m11 * m20 * i[7] + m00 * m21 * i[6] +
              m01 * m20 * i[5] + m00 * m11 * i[4] + m01 * m10 * i[3] +
              m20 * m21 * i[2] + m10 * m11 * i[1] + m00 * m01 * i[0],
          m12 * m20 * i[8] + m10 * m22 * i[7] + m02 * m20 * i[6] +
              m00 * m22 * i[5] + m02 * m10 * i[4] + m00 * m12 * i[3] +
              m20 * m22 * i[2] + m10 * m12 * i[1] + m00 * m02 * i[0],
          m10 * m22 * i[8] + m12 * m20 * i[7] + m00 * m22 * i[6] +
              m02 * m20 * i[5] + m00 * m12 * i[4] + m02 * m10 * i[3] +
              m20 * m22 * i[2] + m10 * m12 * i[1] + m00 * m02 * i[0],
          m12 * m21 * i[8] + m11 * m22 * i[7] + m02 * m21 * i[6] +
              m01 * m22 * i[5] + m02 * m11 * i[4] + m01 * m12 * i[3] +
              m21 * m22 * i[2] + m11 * m12 * i[1] + m01 * m02 * i[0],
          m11 * m22 * i[8] + m12 * m21 * i[7] + m01 * m22 * i[6] +
              m02 * m21 * i[5] + m01 * m12 * i[4] + m02 * m11 * i[3] +
              m21 * m22 * i[2] + m12 * m12 * i[1] + m01 * m02 * i[0],
      };
      o[0] = nv[0];
      o[1] = nv[1];
      o[2] = nv[2];
      o[3] = nv[3];
      o[4] = nv[4];
      o[5] = nv[5];
      o[6] = nv[6];
      o[7] = nv[7];
      o[8] = nv[8];
    }  // end of Rotation3D::rotateTensor

    void changeBasis(real* const v,
                     const std::vector<Variable>& vs,
                     const Hypothesis h,
                     const real* const m) {
      changeBasis(v, v, vs, h, m);
    }  // end of changeBasis

    void changeBasis(real* const o,
                     const real* const i,
                     const std::vector<Variable>& vs,
                     const Hypothesis h,
                     const real* const m) {
      const auto n = getSpaceDimension(h);
      if (n == 1u) {
      } else if (n == 2u) {
        Rotation2D r(m);
        changeBasis(o, i, vs, h, r);
      } else if (n == 3u) {
        Rotation3D r(m);
        changeBasis(o, i, vs, h, r);
      } else {
        raise("changeBasis: unsupported dimension");
      }
    }  // end of changeBasis

  }  // end of namespace behaviour

}  // end of namespace mgis
