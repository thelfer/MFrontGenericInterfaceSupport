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
#include "MGIS/MatrixView.hxx"
#include "MGIS/Behaviour/ChangeBasis.hxx"

namespace mgis {

  namespace behaviour {

    Rotation2D::Rotation2D(const real* const m)
        : m00(m[0]),
          m01(m[1]),
          m10(m[3]),
          m11(m[4]) {}  // end of Rotation2D::Rotation2D

    Rotation2D::Rotation2D(const real v00,
                           const real v01,
                           const real v10,
                           const real v11)
        : m00(v00),
          m01(v01),
          m10(v10),
          m11(v11) {}  // end of Rotation2D::Rotation2D

    Rotation2D Rotation2D::transpose() const {
      return Rotation2D(this->m00, this->m10, this->m01, this->m11);
    }  // end of Rotation2D::transopose

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

    void Rotation2D::buildVectorRotationOperator(mgis::MatrixView& m) const {
      m(0, 0) = this->m00;
      m(0, 1) = this->m01;
      m(1, 0) = this->m10;
      m(1, 1) = this->m11;
    }  // end of Rotation2D::buildVectorRotationOperator

    void Rotation2D::buildStensorRotationOperator(mgis::MatrixView& m) const {
      constexpr const auto cste = Cste::sqrt2;
      constexpr const auto zero = real{0};
      constexpr const auto one = real{1};
      m(0, 0) = this->m00 * this->m00;
      m(0, 1) = this->m10 * this->m10;
      m(0, 2) = zero;
      m(0, 3) = cste * this->m00 * this->m10;
      m(1, 0) = this->m01 * this->m01;
      m(1, 1) = this->m11 * this->m11;
      m(1, 2) = zero;
      m(1, 3) = cste * this->m01 * this->m11;
      m(2, 0) = zero;
      m(2, 1) = zero;
      m(2, 2) = one;
      m(2, 3) = zero;
      m(3, 0) = cste * this->m00 * this->m01;
      m(3, 1) = cste * this->m10 * this->m11;
      m(3, 2) = zero;
      m(3, 3) = this->m00 * this->m11 + this->m01 * this->m10;
    }  // end of Rotation2D::buildStensorRotationOperator

    void Rotation2D::buildTensorRotationOperator(mgis::MatrixView& m) const {
      constexpr const auto zero = real{0};
      constexpr const auto one = real{1};
      m(0, 4) = this->m00 * this->m10;
      m(0, 3) = this->m00 * this->m10;
      m(0, 2) = zero;
      m(0, 1) = this->m10 * this->m10;
      m(0, 0) = this->m00 * this->m00;
      m(1, 4) = this->m01 * this->m11;
      m(1, 3) = this->m01 * this->m11;
      m(1, 2) = zero;
      m(1, 1) = this->m11 * this->m11;
      m(1, 0) = this->m01 * this->m01;
      m(3, 4) = this->m01 * this->m10;
      m(3, 3) = this->m00 * this->m11;
      m(3, 2) = zero;
      m(3, 1) = this->m10 * this->m11;
      m(3, 0) = this->m00 * this->m01;
      m(4, 4) = this->m00 * this->m11;
      m(4, 3) = this->m01 * this->m10;
      m(4, 2) = zero;
      m(4, 1) = this->m10 * this->m11;
      m(4, 0) = this->m00 * this->m01;
      m(2, 4) = zero;
      m(2, 3) = zero;
      m(2, 2) = one;
      m(2, 1) = zero;
      m(2, 0) = zero;
    }  // end of Rotation2D::buildTensorRotationOperator

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

    Rotation3D::Rotation3D(const real v00,
                           const real v01,
                           const real v02,
                           const real v10,
                           const real v11,
                           const real v12,
                           const real v20,
                           const real v21,
                           const real v22)
        : m00(v00),
          m01(v01),
          m02(v02),
          m10(v10),
          m11(v11),
          m12(v12),
          m20(v20),
          m21(v21),
          m22(v22) {}  // end of Rotation3D::Rotation3D

    Rotation3D Rotation3D::transpose() const {
      return Rotation3D(this->m00, this->m10, this->m20, this->m01, this->m11,
                        this->m21, this->m02, this->m12, this->m22);
    }  // end of Rotation3D::transopose

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

    void Rotation3D::buildVectorRotationOperator(mgis::MatrixView& m) const {
      m(0, 0) = this->m00;
      m(0, 1) = this->m01;
      m(0, 2) = this->m02;
      m(1, 0) = this->m10;
      m(1, 1) = this->m11;
      m(1, 2) = this->m12;
      m(2, 0) = this->m20;
      m(2, 1) = this->m21;
      m(2, 2) = this->m22;
    }  // end of Rotation3D::buildVectorRotationOperator

    void Rotation3D::buildStensorRotationOperator(mgis::MatrixView& m) const {
      constexpr const auto cste = Cste::sqrt2;
      m(0, 0) = this->m00 * this->m00;
      m(0, 1) = this->m10 * this->m10;
      m(0, 2) = this->m20 * this->m20;
      m(0, 3) = cste * this->m00 * this->m10;
      m(0, 4) = cste * this->m00 * this->m20;
      m(0, 5) = cste * this->m10 * this->m20;
      m(1, 0) = this->m01 * this->m01;
      m(1, 1) = this->m11 * this->m11;
      m(1, 2) = this->m21 * this->m21;
      m(1, 3) = cste * this->m01 * this->m11;
      m(1, 4) = cste * this->m01 * this->m21;
      m(1, 5) = cste * this->m11 * this->m21;
      m(2, 0) = this->m02 * this->m02;
      m(2, 1) = this->m12 * this->m12;
      m(2, 2) = this->m22 * this->m22;
      m(2, 3) = cste * this->m02 * this->m12;
      m(2, 4) = cste * this->m02 * this->m22;
      m(2, 5) = cste * this->m12 * this->m22;
      m(3, 0) = cste * this->m00 * this->m01;
      m(3, 1) = cste * this->m10 * this->m11;
      m(3, 2) = cste * this->m20 * this->m21;
      m(3, 3) = this->m00 * this->m11 + this->m01 * this->m10;
      m(3, 4) = this->m00 * this->m21 + this->m01 * this->m20;
      m(3, 5) = this->m10 * this->m21 + this->m11 * this->m20;
      m(4, 0) = cste * this->m00 * this->m02;
      m(4, 1) = cste * this->m10 * this->m12;
      m(4, 2) = cste * this->m20 * this->m22;
      m(4, 3) = this->m00 * this->m12 + this->m02 * this->m10;
      m(4, 4) = this->m00 * this->m22 + this->m02 * this->m20;
      m(4, 5) = this->m10 * this->m22 + this->m12 * this->m20;
      m(5, 0) = cste * this->m01 * this->m02;
      m(5, 1) = cste * this->m11 * this->m12;
      m(5, 2) = cste * this->m21 * this->m22;
      m(5, 3) = this->m01 * this->m12 + this->m02 * this->m11;
      m(5, 4) = this->m01 * this->m22 + this->m02 * this->m21;
      m(5, 5) = this->m11 * this->m22 + this->m12 * this->m21;
    }  // end of Rotation3D::buildStensorRotationOperator

    void Rotation3D::buildTensorRotationOperator(mgis::MatrixView& m) const {
      m(0, 0) = this->m00 * this->m00;
      m(0, 1) = this->m10 * this->m10;
      m(0, 2) = this->m20 * this->m20;
      m(0, 3) = this->m00 * this->m10;
      m(0, 4) = this->m00 * this->m10;
      m(0, 5) = this->m00 * this->m20;
      m(0, 6) = this->m00 * this->m20;
      m(0, 7) = this->m10 * this->m20;
      m(0, 8) = this->m10 * this->m20;
      m(1, 0) = this->m01 * this->m01;
      m(1, 1) = this->m11 * this->m11;
      m(1, 2) = this->m21 * this->m21;
      m(1, 3) = this->m01 * this->m11;
      m(1, 4) = this->m01 * this->m11;
      m(1, 5) = this->m01 * this->m21;
      m(1, 6) = this->m01 * this->m21;
      m(1, 7) = this->m11 * this->m21;
      m(1, 8) = this->m11 * this->m21;
      m(2, 0) = this->m02 * this->m02;
      m(2, 1) = this->m12 * this->m12;
      m(2, 2) = this->m22 * this->m22;
      m(2, 3) = this->m02 * this->m12;
      m(2, 4) = this->m02 * this->m12;
      m(2, 5) = this->m02 * this->m22;
      m(2, 6) = this->m02 * this->m22;
      m(2, 7) = this->m12 * this->m22;
      m(2, 8) = this->m12 * this->m22;
      m(3, 8) = this->m11 * this->m20;
      m(3, 7) = this->m10 * this->m21;
      m(3, 6) = this->m01 * this->m20;
      m(3, 5) = this->m00 * this->m21;
      m(3, 4) = this->m01 * this->m10;
      m(3, 3) = this->m00 * this->m11;
      m(3, 2) = this->m20 * this->m21;
      m(3, 1) = this->m10 * this->m11;
      m(3, 0) = this->m00 * this->m01;
      m(4, 8) = this->m10 * this->m21;
      m(4, 7) = this->m11 * this->m20;
      m(4, 6) = this->m00 * this->m21;
      m(4, 5) = this->m01 * this->m20;
      m(4, 4) = this->m00 * this->m11;
      m(4, 3) = this->m01 * this->m10;
      m(4, 2) = this->m20 * this->m21;
      m(4, 1) = this->m10 * this->m11;
      m(4, 0) = this->m00 * this->m01;
      m(5, 8) = this->m12 * this->m20;
      m(5, 7) = this->m10 * this->m22;
      m(5, 6) = this->m02 * this->m20;
      m(5, 5) = this->m00 * this->m22;
      m(5, 4) = this->m02 * this->m10;
      m(5, 3) = this->m00 * this->m12;
      m(5, 2) = this->m20 * this->m22;
      m(5, 1) = this->m10 * this->m12;
      m(5, 0) = this->m00 * this->m02;
      m(6, 8) = this->m10 * this->m22;
      m(6, 7) = this->m12 * this->m20;
      m(6, 6) = this->m00 * this->m22;
      m(6, 5) = this->m02 * this->m20;
      m(6, 4) = this->m00 * this->m12;
      m(6, 3) = this->m02 * this->m10;
      m(6, 2) = this->m20 * this->m22;
      m(6, 1) = this->m10 * this->m12;
      m(6, 0) = this->m00 * this->m02;
      m(7, 8) = this->m12 * this->m21;
      m(7, 7) = this->m11 * this->m22;
      m(7, 6) = this->m02 * this->m21;
      m(7, 5) = this->m01 * this->m22;
      m(7, 4) = this->m02 * this->m11;
      m(7, 3) = this->m01 * this->m12;
      m(7, 2) = this->m21 * this->m22;
      m(7, 1) = this->m11 * this->m12;
      m(7, 0) = this->m01 * this->m02;
      m(8, 8) = this->m11 * this->m22;
      m(8, 7) = this->m12 * this->m21;
      m(8, 6) = this->m01 * this->m22;
      m(8, 5) = this->m02 * this->m21;
      m(8, 4) = this->m01 * this->m12;
      m(8, 3) = this->m02 * this->m11;
      m(8, 2) = this->m21 * this->m22;
      m(8, 1) = this->m11 * this->m12;
      m(8, 0) = this->m01 * this->m02;
    }  // end of Rotation3D::buildTensorRotationOperator

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
