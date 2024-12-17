/*!
 * \file   src/FiniteStrainSupport.cxx
 * \brief
 * \author Thomas Helfer
 * \date   25/01/2019
 * \copyright (C) Copyright Thomas Helfer 2018-2019.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#include "MGIS/Raise.hxx"
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/BehaviourData.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/FiniteStrainSupport.hxx"

namespace mgis::behaviour {

  void convertFiniteStrainStress_PK1_2D(real* const P,
                                        const real* const F,
                                        const real* const s) {
    constexpr const real cste = 1.41421356237309504880;
    constexpr const real icste = 0.70710678118654752440;
    P[0] = -(s[3] * F[2] * F[3] - cste * s[0] * F[1] * F[2]) * icste;
    P[1] = -(s[3] * F[2] * F[4] - cste * F[0] * s[1] * F[2]) * icste;
    P[2] = F[0] * s[2] * F[1] - s[2] * F[3] * F[4];
    P[3] = -(cste * s[0] * F[2] * F[4] - F[0] * s[3] * F[2]) * icste;
    P[4] = -(cste * s[1] * F[2] * F[3] - s[3] * F[1] * F[2]) * icste;
  }  // end of convertFiniteStrainStress_PK1_2D

  void convertFiniteStrainStress_PK1_3D(real* const P,
                                        const real* const F,
                                        const real* const s) {
    constexpr const real cste = 1.41421356237309504880;
    P[0] = -((2 * s[0] * F[7] - cste * s[3] * F[5]) * F[8] -
             cste * s[4] * F[3] * F[7] + cste * s[4] * F[1] * F[5] +
             cste * s[3] * F[2] * F[3] - 2 * s[0] * F[1] * F[2]) /
           2;
    P[1] = ((cste * s[3] * F[6] - cste * F[0] * s[5]) * F[7] -
            2 * s[1] * F[5] * F[6] + cste * s[5] * F[4] * F[5] -
            cste * s[3] * F[2] * F[4] + 2 * F[0] * s[1] * F[2]) /
           2;
    P[2] = ((cste * s[4] * F[4] - cste * F[0] * s[5]) * F[8] +
            (cste * s[5] * F[3] - cste * s[4] * F[1]) * F[6] -
            2 * s[2] * F[3] * F[4] + 2 * F[0] * s[2] * F[1]) /
           2;
    P[3] = ((2 * s[0] * F[6] - cste * F[0] * s[4]) * F[7] -
            cste * s[3] * F[5] * F[6] + cste * s[4] * F[4] * F[5] -
            2 * s[0] * F[2] * F[4] + cste * F[0] * s[3] * F[2]) /
           2;
    P[4] = -((cste * s[3] * F[7] - 2 * s[1] * F[5]) * F[8] -
             cste * s[5] * F[3] * F[7] + cste * s[5] * F[1] * F[5] +
             2 * s[1] * F[2] * F[3] - cste * s[3] * F[1] * F[2]) /
           2;
    P[5] = ((2 * s[0] * F[4] - cste * F[0] * s[3]) * F[8] +
            (cste * s[3] * F[3] - 2 * s[0] * F[1]) * F[6] -
            cste * s[4] * F[3] * F[4] + cste * F[0] * s[4] * F[1]) /
           2;
    P[6] = -((cste * s[4] * F[7] - cste * s[5] * F[5]) * F[8] -
             2 * s[2] * F[3] * F[7] + 2 * s[2] * F[1] * F[5] +
             cste * s[5] * F[2] * F[3] - cste * s[4] * F[1] * F[2]) /
           2;
    P[7] = ((cste * s[3] * F[4] - 2 * F[0] * s[1]) * F[8] +
            (2 * s[1] * F[3] - cste * s[3] * F[1]) * F[6] -
            cste * s[5] * F[3] * F[4] + cste * F[0] * s[5] * F[1]) /
           2;
    P[8] = ((cste * s[4] * F[6] - 2 * F[0] * s[2]) * F[7] -
            cste * s[5] * F[5] * F[6] + 2 * s[2] * F[4] * F[5] -
            cste * s[4] * F[2] * F[4] + cste * F[0] * s[5] * F[2]) /
           2;
  }  // end of convertFiniteStrainStress_PK1_3D

  static void convertFiniteStrainStress_PK1_2D(std::span<real>& Ps,
                                               const MaterialDataManager& m,
                                               const mgis::size_type b,
                                               const mgis::size_type e) {
    auto* const P = Ps.data();
    const auto* const F = m.s1.gradients.data();
    const auto* const s = m.s1.thermodynamic_forces.data();
    // symmetric tensor size
    const auto ss = getStensorSize(Hypothesis::PLANESTRAIN);
    // non symmetric tensor size
    const auto ts = getTensorSize(Hypothesis::PLANESTRAIN);
    for (auto i = b; i != e; ++i) {
      auto* const P_l = P + ts * i;
      const auto* const F_l = F + ts * i;
      const auto* const s_l = s + ss * i;
      convertFiniteStrainStress_PK1_2D(P_l, F_l, s_l);
    }
  }  // end of convertFiniteStrainStress_PK1_2D

  static void convertFiniteStrainStress_PK1_3D(std::span<real>& Ps,
                                               const MaterialDataManager& m,
                                               const mgis::size_type b,
                                               const mgis::size_type e) {
    auto* const P = Ps.data();
    const auto* const F = m.s1.gradients.data();
    const auto* const s = m.s1.thermodynamic_forces.data();
    // symmetric tensor size
    const auto ss = getStensorSize(Hypothesis::TRIDIMENSIONAL);
    // non symmetric tensor size
    const auto ts = getTensorSize(Hypothesis::TRIDIMENSIONAL);
    for (auto i = b; i != e; ++i) {
      auto* const P_l = P + ts * i;
      const auto* const F_l = F + ts * i;
      const auto* const s_l = s + ss * i;
      convertFiniteStrainStress_PK1_3D(P_l, F_l, s_l);
    }
  }  // end of convertFiniteStrainStress_PK1_3D

  static void convertFiniteStrainStress_PK1_2D(std::span<real>& s,
                                               const MaterialDataManager& m) {
    // check behaviour type
    if (m.b.btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      mgis::raise(
          "convertFiniteStrainStress: "
          "unsupported behaviour type");
    }
    // non symmetric tensor size
    const auto ts = getTensorSize(m.b.hypothesis);
    // check K size
    if (s.size() != m.n * ts) {
      mgis::raise(
          "convertFiniteStrainStress: "
          "unsupported tangent operator");
    }
    convertFiniteStrainStress_PK1_2D(s, m, 0, m.n);
  }  // end of convertFiniteStrainStress_PK1_2D

  static void convertFiniteStrainStress_PK1_3D(std::span<real>& s,
                                               const MaterialDataManager& m) {
    // check behaviour type
    if (m.b.btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      mgis::raise(
          "convertFiniteStrainStress: "
          "unsupported behaviour type");
    }
    // non symmetric tensor size
    const auto ts = getTensorSize(m.b.hypothesis);
    // check K size
    if (s.size() != m.n * ts) {
      mgis::raise(
          "convertFiniteStrainStress: "
          "unsupported tangent operator");
    }
    convertFiniteStrainStress_PK1_3D(s, m, 0, m.n);
  }  // end of convertFiniteStrainStress_PK1_3D

  void convertFiniteStrainStress(std::span<real>& s,
                                 const MaterialDataManager& m,
                                 const FiniteStrainStress t) {
    const auto h = m.b.hypothesis;
    if (t == FiniteStrainStress::PK1) {
      if (h == Hypothesis::TRIDIMENSIONAL) {
        convertFiniteStrainStress_PK1_3D(s, m);
      } else if ((h == Hypothesis::AXISYMMETRICAL) ||
                 (h == Hypothesis::PLANESTRAIN) ||
                 (h == Hypothesis::GENERALISEDPLANESTRAIN)) {
        convertFiniteStrainStress_PK1_2D(s, m);
      } else {
        mgis::raise(
            "convertFiniteStrainStress: "
            "unsupported hypothesis");
      }
    } else {
      mgis::raise(
          "convertFiniteStrainStress: "
          "unsupported tangent operator");
    }
  }  // end of convertFiniteStrainStress

  static void convertFiniteStrainStress_PK1_2D(std::span<real>& P,
                                               const BehaviourData& d) {
    // check behaviour type
    if (d.s1.b.btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      mgis::raise(
          "convertFiniteStrainStress: "
          "unsupported behaviour type");
    }
    // non symmetric tensor size
    const auto ts = getTensorSize(d.s1.b.hypothesis);
    // check K size
    if (P.size() != ts) {
      mgis::raise(
          "convertFiniteStrainStress: "
          "unsupported tangent operator");
    }
    convertFiniteStrainStress_PK1_2D(P.data(), d.s1.gradients.data(),
                                     d.s1.thermodynamic_forces.data());
  }  // end of convertFiniteStrainStress_PK1_2D

  static void convertFiniteStrainStress_PK1_3D(std::span<real>& P,
                                               const BehaviourData& d) {
    // check behaviour type
    if (d.s1.b.btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      mgis::raise(
          "convertFiniteStrainStress: "
          "unsupported behaviour type");
    }
    // non symmetric tensor size
    const auto ts = getTensorSize(d.s1.b.hypothesis);
    // check K size
    if (P.size() != ts) {
      mgis::raise(
          "convertFiniteStrainStress: "
          "unsupported tangent operator");
    }
    convertFiniteStrainStress_PK1_3D(P.data(), d.s1.gradients.data(),
                                     d.s1.thermodynamic_forces.data());
  }  // end of convertFiniteStrainStress_PK1_3D

  void convertFiniteStrainStress(std::span<real>& s,
                                 const BehaviourData& d,
                                 const FiniteStrainStress t) {
    const auto h = d.s1.b.hypothesis;
    if (t == FiniteStrainStress::PK1) {
      if (h == Hypothesis::TRIDIMENSIONAL) {
        convertFiniteStrainStress_PK1_3D(s, d);
      } else if ((h == Hypothesis::AXISYMMETRICAL) ||
                 (h == Hypothesis::PLANESTRAIN) ||
                 (h == Hypothesis::GENERALISEDPLANESTRAIN)) {
        convertFiniteStrainStress_PK1_2D(s, d);
      } else {
        mgis::raise(
            "convertFiniteStrainStress: "
            "unsupported hypothesis");
      }
    } else {
      mgis::raise(
          "convertFiniteStrainStress: "
          "unsupported tangent operator");
    }
  }  // end of convertFiniteStrainStress

  static void convertFiniteStrainTangentOperator_PK1_2D(mgis::real* const dP,
                                                        const real* const ds,
                                                        const real* const F,
                                                        const real* const s) {
    constexpr const real cste = 1.41421356237309504880;
    constexpr const real icste = 0.70710678118654752440;
    // diff(P[1],F[0])
    dP[0] = -(ds[15] * F[2] * F[3] - cste * ds[0] * F[1] * F[2]) * icste;
    dP[5] = -(ds[15] * F[2] * F[4] - cste * F[0] * ds[5] * F[2] -
              cste * s[1] * F[2]) *
            icste;
    dP[10] = (-ds[10] * F[3] * F[4]) + F[0] * ds[10] * F[1] + s[2] * F[1];
    dP[15] =
        -(cste * ds[0] * F[2] * F[4] - F[0] * ds[15] * F[2] - s[3] * F[2]) *
        icste;
    dP[20] = -(cste * ds[5] * F[2] * F[3] - ds[15] * F[1] * F[2]) * icste;
    // diff(P[1],F[1])
    dP[1] = -(ds[16] * F[2] * F[3] - cste * ds[1] * F[1] * F[2] -
              cste * s[0] * F[2]) *
            icste;
    dP[6] = -(ds[16] * F[2] * F[4] - cste * F[0] * ds[6] * F[2]) * icste;
    dP[11] = (-ds[11] * F[3] * F[4]) + F[0] * ds[11] * F[1] + F[0] * s[2];
    dP[16] = -(cste * ds[1] * F[2] * F[4] - F[0] * ds[16] * F[2]) * icste;
    dP[21] =
        -(cste * ds[6] * F[2] * F[3] - ds[16] * F[1] * F[2] - s[3] * F[2]) *
        icste;
    // diff(P[1],F[2])
    dP[2] = -(ds[17] * F[2] * F[3] + s[3] * F[3] - cste * ds[2] * F[1] * F[2] -
              cste * s[0] * F[1]) *
            icste;
    dP[7] = -(ds[17] * F[2] * F[4] + s[3] * F[4] - cste * F[0] * ds[7] * F[2] -
              cste * F[0] * s[1]) *
            icste;
    dP[12] = F[0] * ds[12] * F[1] - ds[12] * F[3] * F[4];
    dP[17] = -(cste * ds[2] * F[2] * F[4] + cste * s[0] * F[4] -
               F[0] * ds[17] * F[2] - F[0] * s[3]) *
             icste;
    dP[22] = -(cste * ds[7] * F[2] * F[3] + cste * s[1] * F[3] -
               ds[17] * F[1] * F[2] - s[3] * F[1]) *
             icste;
    // diff(P[1],F[3])
    dP[3] = -(ds[18] * F[2] * F[3] - cste * ds[3] * F[1] * F[2] + s[3] * F[2]) *
            icste;
    dP[8] = -(ds[18] * F[2] * F[4] - cste * F[0] * ds[8] * F[2]) * icste;
    dP[13] = (-ds[13] * F[3] * F[4]) - s[2] * F[4] + F[0] * ds[13] * F[1];
    dP[18] = -(cste * ds[3] * F[2] * F[4] - F[0] * ds[18] * F[2]) * icste;
    dP[23] = -(cste * ds[8] * F[2] * F[3] - ds[18] * F[1] * F[2] +
               cste * s[1] * F[2]) *
             icste;
    // f90(diff(P[1],F[4]));
    dP[4] = -(ds[19] * F[2] * F[3] - cste * ds[4] * F[1] * F[2]) * icste;
    dP[9] = -(ds[19] * F[2] * F[4] - cste * F[0] * ds[9] * F[2] + s[3] * F[2]) *
            icste;
    dP[14] = (-ds[14] * F[3] * F[4]) - s[2] * F[3] + F[0] * ds[14] * F[1];
    dP[19] = -(cste * ds[4] * F[2] * F[4] - F[0] * ds[19] * F[2] +
               cste * s[0] * F[2]) *
             icste;
    dP[24] = -(cste * ds[9] * F[2] * F[3] - ds[19] * F[1] * F[2]) * icste;
  }  // end of convertFiniteStrainTangentOperator_PK1_2D

  static void convertFiniteStrainTangentOperator_PK1_3D(mgis::real* const dP,
                                                        const real* const ds,
                                                        const real* const F,
                                                        const real* const s) {
    constexpr const real cste = 1.41421356237309504880;
    //(%i15) f90(diff(P[1],F_0));
    dP[0] = -((2 * ds[0] * F[7] - cste * ds[27] * F[5]) * F[8] -
              cste * ds[36] * F[3] * F[7] + cste * ds[36] * F[1] * F[5] +
              cste * ds[27] * F[2] * F[3] - 2 * ds[0] * F[1] * F[2]) /
            2;
    dP[9] =
        ((cste * ds[27] * F[6] - cste * F[0] * ds[45] - cste * s[5]) * F[7] -
         2 * ds[9] * F[5] * F[6] + cste * ds[45] * F[4] * F[5] -
         cste * ds[27] * F[2] * F[4] + 2 * F[0] * ds[9] * F[2] +
         2 * s[1] * F[2]) /
        2;
    dP[18] =
        ((cste * ds[36] * F[4] - cste * F[0] * ds[45] - cste * s[5]) * F[8] +
         (cste * ds[45] * F[3] - cste * ds[36] * F[1]) * F[6] -
         2 * ds[18] * F[3] * F[4] + 2 * F[0] * ds[18] * F[1] +
         2 * s[2] * F[1]) /
        2;
    dP[27] = ((2 * ds[0] * F[6] - cste * F[0] * ds[36] - cste * s[4]) * F[7] -
              cste * ds[27] * F[5] * F[6] + cste * ds[36] * F[4] * F[5] -
              2 * ds[0] * F[2] * F[4] + cste * F[0] * ds[27] * F[2] +
              cste * s[3] * F[2]) /
             2;
    dP[36] = -((cste * ds[27] * F[7] - 2 * ds[9] * F[5]) * F[8] -
               cste * ds[45] * F[3] * F[7] + cste * ds[45] * F[1] * F[5] +
               2 * ds[9] * F[2] * F[3] - cste * ds[27] * F[1] * F[2]) /
             2;
    dP[45] = ((2 * ds[0] * F[4] - cste * F[0] * ds[27] - cste * s[3]) * F[8] +
              (cste * ds[27] * F[3] - 2 * ds[0] * F[1]) * F[6] -
              cste * ds[36] * F[3] * F[4] + cste * F[0] * ds[36] * F[1] +
              cste * s[4] * F[1]) /
             2;
    dP[54] = -((cste * ds[36] * F[7] - cste * ds[45] * F[5]) * F[8] -
               2 * ds[18] * F[3] * F[7] + 2 * ds[18] * F[1] * F[5] +
               cste * ds[45] * F[2] * F[3] - cste * ds[36] * F[1] * F[2]) /
             2;
    dP[63] = ((cste * ds[27] * F[4] - 2 * F[0] * ds[9] - 2 * s[1]) * F[8] +
              (2 * ds[9] * F[3] - cste * ds[27] * F[1]) * F[6] -
              cste * ds[45] * F[3] * F[4] + cste * F[0] * ds[45] * F[1] +
              cste * s[5] * F[1]) /
             2;
    dP[72] = ((cste * ds[36] * F[6] - 2 * F[0] * ds[18] - 2 * s[2]) * F[7] -
              cste * ds[45] * F[5] * F[6] + 2 * ds[18] * F[4] * F[5] -
              cste * ds[36] * F[2] * F[4] + cste * F[0] * ds[45] * F[2] +
              cste * s[5] * F[2]) /
             2;
    // (%i16) f90(diff(P[1],F[1]));
    dP[1] = -((2 * ds[1] * F[7] - cste * ds[28] * F[5]) * F[8] -
              cste * ds[37] * F[3] * F[7] + cste * ds[37] * F[1] * F[5] +
              cste * s[4] * F[5] + cste * ds[28] * F[2] * F[3] -
              2 * ds[1] * F[1] * F[2] - 2 * s[0] * F[2]) /
            2;
    dP[10] = ((cste * ds[28] * F[6] - cste * F[0] * ds[46]) * F[7] -
              2 * ds[10] * F[5] * F[6] + cste * ds[46] * F[4] * F[5] -
              cste * ds[28] * F[2] * F[4] + 2 * F[0] * ds[10] * F[2]) /
             2;
    dP[19] =
        ((cste * ds[37] * F[4] - cste * F[0] * ds[46]) * F[8] +
         (cste * ds[46] * F[3] - cste * ds[37] * F[1] - cste * s[4]) * F[6] -
         2 * ds[19] * F[3] * F[4] + 2 * F[0] * ds[19] * F[1] +
         2 * F[0] * s[2]) /
        2;
    dP[28] = ((2 * ds[1] * F[6] - cste * F[0] * ds[37]) * F[7] -
              cste * ds[28] * F[5] * F[6] + cste * ds[37] * F[4] * F[5] -
              2 * ds[1] * F[2] * F[4] + cste * F[0] * ds[28] * F[2]) /
             2;
    dP[37] = -((cste * ds[28] * F[7] - 2 * ds[10] * F[5]) * F[8] -
               cste * ds[46] * F[3] * F[7] + cste * ds[46] * F[1] * F[5] +
               cste * s[5] * F[5] + 2 * ds[10] * F[2] * F[3] -
               cste * ds[28] * F[1] * F[2] - cste * s[3] * F[2]) /
             2;
    dP[46] = ((2 * ds[1] * F[4] - cste * F[0] * ds[28]) * F[8] +
              (cste * ds[28] * F[3] - 2 * ds[1] * F[1] - 2 * s[0]) * F[6] -
              cste * ds[37] * F[3] * F[4] + cste * F[0] * ds[37] * F[1] +
              cste * F[0] * s[4]) /
             2;
    dP[55] = -((cste * ds[37] * F[7] - cste * ds[46] * F[5]) * F[8] -
               2 * ds[19] * F[3] * F[7] + 2 * ds[19] * F[1] * F[5] +
               2 * s[2] * F[5] + cste * ds[46] * F[2] * F[3] -
               cste * ds[37] * F[1] * F[2] - cste * s[4] * F[2]) /
             2;
    dP[64] = ((cste * ds[28] * F[4] - 2 * F[0] * ds[10]) * F[8] +
              (2 * ds[10] * F[3] - cste * ds[28] * F[1] - cste * s[3]) * F[6] -
              cste * ds[46] * F[3] * F[4] + cste * F[0] * ds[46] * F[1] +
              cste * F[0] * s[5]) /
             2;
    dP[73] = ((cste * ds[37] * F[6] - 2 * F[0] * ds[19]) * F[7] -
              cste * ds[46] * F[5] * F[6] + 2 * ds[19] * F[4] * F[5] -
              cste * ds[37] * F[2] * F[4] + cste * F[0] * ds[46] * F[2]) /
             2;
    // (%i17) f90(diff(P[1],F[2]));
    dP[2] = -((2 * ds[2] * F[7] - cste * ds[29] * F[5]) * F[8] -
              cste * ds[38] * F[3] * F[7] + cste * ds[38] * F[1] * F[5] +
              cste * ds[29] * F[2] * F[3] + cste * s[3] * F[3] -
              2 * ds[2] * F[1] * F[2] - 2 * s[0] * F[1]) /
            2;
    dP[11] = ((cste * ds[29] * F[6] - cste * F[0] * ds[47]) * F[7] -
              2 * ds[11] * F[5] * F[6] + cste * ds[47] * F[4] * F[5] -
              cste * ds[29] * F[2] * F[4] - cste * s[3] * F[4] +
              2 * F[0] * ds[11] * F[2] + 2 * F[0] * s[1]) /
             2;
    dP[20] = ((cste * ds[38] * F[4] - cste * F[0] * ds[47]) * F[8] +
              (cste * ds[47] * F[3] - cste * ds[38] * F[1]) * F[6] -
              2 * ds[20] * F[3] * F[4] + 2 * F[0] * ds[20] * F[1]) /
             2;
    dP[29] = ((2 * ds[2] * F[6] - cste * F[0] * ds[38]) * F[7] -
              cste * ds[29] * F[5] * F[6] + cste * ds[38] * F[4] * F[5] -
              2 * ds[2] * F[2] * F[4] - 2 * s[0] * F[4] +
              cste * F[0] * ds[29] * F[2] + cste * F[0] * s[3]) /
             2;
    dP[38] = -((cste * ds[29] * F[7] - 2 * ds[11] * F[5]) * F[8] -
               cste * ds[47] * F[3] * F[7] + cste * ds[47] * F[1] * F[5] +
               2 * ds[11] * F[2] * F[3] + 2 * s[1] * F[3] -
               cste * ds[29] * F[1] * F[2] - cste * s[3] * F[1]) /
             2;
    dP[47] = ((2 * ds[2] * F[4] - cste * F[0] * ds[29]) * F[8] +
              (cste * ds[29] * F[3] - 2 * ds[2] * F[1]) * F[6] -
              cste * ds[38] * F[3] * F[4] + cste * F[0] * ds[38] * F[1]) /
             2;
    dP[56] = -((cste * ds[38] * F[7] - cste * ds[47] * F[5]) * F[8] -
               2 * ds[20] * F[3] * F[7] + 2 * ds[20] * F[1] * F[5] +
               cste * ds[47] * F[2] * F[3] + cste * s[5] * F[3] -
               cste * ds[38] * F[1] * F[2] - cste * s[4] * F[1]) /
             2;
    dP[65] = ((cste * ds[29] * F[4] - 2 * F[0] * ds[11]) * F[8] +
              (2 * ds[11] * F[3] - cste * ds[29] * F[1]) * F[6] -
              cste * ds[47] * F[3] * F[4] + cste * F[0] * ds[47] * F[1]) /
             2;
    dP[74] = ((cste * ds[38] * F[6] - 2 * F[0] * ds[20]) * F[7] -
              cste * ds[47] * F[5] * F[6] + 2 * ds[20] * F[4] * F[5] -
              cste * ds[38] * F[2] * F[4] - cste * s[4] * F[4] +
              cste * F[0] * ds[47] * F[2] + cste * F[0] * s[5]) /
             2;
    // (%i18) f90(diff(P[1],F[3]));
    dP[3] = -((2 * ds[3] * F[7] - cste * ds[30] * F[5]) * F[8] -
              cste * ds[39] * F[3] * F[7] - cste * s[4] * F[7] +
              cste * ds[39] * F[1] * F[5] + cste * ds[30] * F[2] * F[3] -
              2 * ds[3] * F[1] * F[2] + cste * s[3] * F[2]) /
            2;
    dP[12] = ((cste * ds[30] * F[6] - cste * F[0] * ds[48]) * F[7] -
              2 * ds[12] * F[5] * F[6] + cste * ds[48] * F[4] * F[5] -
              cste * ds[30] * F[2] * F[4] + 2 * F[0] * ds[12] * F[2]) /
             2;
    dP[21] =
        ((cste * ds[39] * F[4] - cste * F[0] * ds[48]) * F[8] +
         (cste * ds[48] * F[3] - cste * ds[39] * F[1] + cste * s[5]) * F[6] -
         2 * ds[21] * F[3] * F[4] - 2 * s[2] * F[4] +
         2 * F[0] * ds[21] * F[1]) /
        2;
    dP[30] = ((2 * ds[3] * F[6] - cste * F[0] * ds[39]) * F[7] -
              cste * ds[30] * F[5] * F[6] + cste * ds[39] * F[4] * F[5] -
              2 * ds[3] * F[2] * F[4] + cste * F[0] * ds[30] * F[2]) /
             2;
    dP[39] = -((cste * ds[30] * F[7] - 2 * ds[12] * F[5]) * F[8] -
               cste * ds[48] * F[3] * F[7] - cste * s[5] * F[7] +
               cste * ds[48] * F[1] * F[5] + 2 * ds[12] * F[2] * F[3] -
               cste * ds[30] * F[1] * F[2] + 2 * s[1] * F[2]) /
             2;
    dP[48] = ((2 * ds[3] * F[4] - cste * F[0] * ds[30]) * F[8] +
              (cste * ds[30] * F[3] - 2 * ds[3] * F[1] + cste * s[3]) * F[6] -
              cste * ds[39] * F[3] * F[4] - cste * s[4] * F[4] +
              cste * F[0] * ds[39] * F[1]) /
             2;
    dP[57] = -((cste * ds[39] * F[7] - cste * ds[48] * F[5]) * F[8] -
               2 * ds[21] * F[3] * F[7] - 2 * s[2] * F[7] +
               2 * ds[21] * F[1] * F[5] + cste * ds[48] * F[2] * F[3] -
               cste * ds[39] * F[1] * F[2] + cste * s[5] * F[2]) /
             2;
    dP[66] = ((cste * ds[30] * F[4] - 2 * F[0] * ds[12]) * F[8] +
              (2 * ds[12] * F[3] - cste * ds[30] * F[1] + 2 * s[1]) * F[6] -
              cste * ds[48] * F[3] * F[4] - cste * s[5] * F[4] +
              cste * F[0] * ds[48] * F[1]) /
             2;
    dP[75] = ((cste * ds[39] * F[6] - 2 * F[0] * ds[21]) * F[7] -
              cste * ds[48] * F[5] * F[6] + 2 * ds[21] * F[4] * F[5] -
              cste * ds[39] * F[2] * F[4] + cste * F[0] * ds[48] * F[2]) /
             2;
    // (%i19) f90(diff(P[1],F[4]));
    dP[4] = -((2 * ds[4] * F[7] - cste * ds[31] * F[5]) * F[8] -
              cste * ds[40] * F[3] * F[7] + cste * ds[40] * F[1] * F[5] +
              cste * ds[31] * F[2] * F[3] - 2 * ds[4] * F[1] * F[2]) /
            2;
    dP[13] = ((cste * ds[31] * F[6] - cste * F[0] * ds[49]) * F[7] -
              2 * ds[13] * F[5] * F[6] + cste * ds[49] * F[4] * F[5] +
              cste * s[5] * F[5] - cste * ds[31] * F[2] * F[4] +
              2 * F[0] * ds[13] * F[2] - cste * s[3] * F[2]) /
             2;
    dP[22] =
        ((cste * ds[40] * F[4] - cste * F[0] * ds[49] + cste * s[4]) * F[8] +
         (cste * ds[49] * F[3] - cste * ds[40] * F[1]) * F[6] -
         2 * ds[22] * F[3] * F[4] - 2 * s[2] * F[3] +
         2 * F[0] * ds[22] * F[1]) /
        2;
    dP[31] = ((2 * ds[4] * F[6] - cste * F[0] * ds[40]) * F[7] -
              cste * ds[31] * F[5] * F[6] + cste * ds[40] * F[4] * F[5] +
              cste * s[4] * F[5] - 2 * ds[4] * F[2] * F[4] +
              cste * F[0] * ds[31] * F[2] - 2 * s[0] * F[2]) /
             2;
    dP[40] = -((cste * ds[31] * F[7] - 2 * ds[13] * F[5]) * F[8] -
               cste * ds[49] * F[3] * F[7] + cste * ds[49] * F[1] * F[5] +
               2 * ds[13] * F[2] * F[3] - cste * ds[31] * F[1] * F[2]) /
             2;
    dP[49] = ((2 * ds[4] * F[4] - cste * F[0] * ds[31] + 2 * s[0]) * F[8] +
              (cste * ds[31] * F[3] - 2 * ds[4] * F[1]) * F[6] -
              cste * ds[40] * F[3] * F[4] - cste * s[4] * F[3] +
              cste * F[0] * ds[40] * F[1]) /
             2;
    dP[58] = -((cste * ds[40] * F[7] - cste * ds[49] * F[5]) * F[8] -
               2 * ds[22] * F[3] * F[7] + 2 * ds[22] * F[1] * F[5] +
               cste * ds[49] * F[2] * F[3] - cste * ds[40] * F[1] * F[2]) /
             2;
    dP[67] = ((cste * ds[31] * F[4] - 2 * F[0] * ds[13] + cste * s[3]) * F[8] +
              (2 * ds[13] * F[3] - cste * ds[31] * F[1]) * F[6] -
              cste * ds[49] * F[3] * F[4] - cste * s[5] * F[3] +
              cste * F[0] * ds[49] * F[1]) /
             2;
    dP[76] = ((cste * ds[40] * F[6] - 2 * F[0] * ds[22]) * F[7] -
              cste * ds[49] * F[5] * F[6] + 2 * ds[22] * F[4] * F[5] +
              2 * s[2] * F[5] - cste * ds[40] * F[2] * F[4] +
              cste * F[0] * ds[49] * F[2] - cste * s[4] * F[2]) /
             2;
    // (%i20) f90(diff(P[1],F[5]));
    dP[5] = -((2 * ds[5] * F[7] - cste * ds[32] * F[5] - cste * s[3]) * F[8] -
              cste * ds[41] * F[3] * F[7] + cste * ds[41] * F[1] * F[5] +
              cste * ds[32] * F[2] * F[3] - 2 * ds[5] * F[1] * F[2] +
              cste * s[4] * F[1]) /
            2;
    dP[14] = ((cste * ds[32] * F[6] - cste * F[0] * ds[50]) * F[7] -
              2 * ds[14] * F[5] * F[6] - 2 * s[1] * F[6] +
              cste * ds[50] * F[4] * F[5] - cste * ds[32] * F[2] * F[4] +
              cste * s[5] * F[4] + 2 * F[0] * ds[14] * F[2]) /
             2;
    dP[23] = ((cste * ds[41] * F[4] - cste * F[0] * ds[50]) * F[8] +
              (cste * ds[50] * F[3] - cste * ds[41] * F[1]) * F[6] -
              2 * ds[23] * F[3] * F[4] + 2 * F[0] * ds[23] * F[1]) /
             2;
    dP[32] = ((2 * ds[5] * F[6] - cste * F[0] * ds[41]) * F[7] -
              cste * ds[32] * F[5] * F[6] - cste * s[3] * F[6] +
              cste * ds[41] * F[4] * F[5] - 2 * ds[5] * F[2] * F[4] +
              cste * s[4] * F[4] + cste * F[0] * ds[32] * F[2]) /
             2;
    dP[41] = -((cste * ds[32] * F[7] - 2 * ds[14] * F[5] - 2 * s[1]) * F[8] -
               cste * ds[50] * F[3] * F[7] + cste * ds[50] * F[1] * F[5] +
               2 * ds[14] * F[2] * F[3] - cste * ds[32] * F[1] * F[2] +
               cste * s[5] * F[1]) /
             2;
    dP[50] = ((2 * ds[5] * F[4] - cste * F[0] * ds[32]) * F[8] +
              (cste * ds[32] * F[3] - 2 * ds[5] * F[1]) * F[6] -
              cste * ds[41] * F[3] * F[4] + cste * F[0] * ds[41] * F[1]) /
             2;
    dP[59] =
        -((cste * ds[41] * F[7] - cste * ds[50] * F[5] - cste * s[5]) * F[8] -
          2 * ds[23] * F[3] * F[7] + 2 * ds[23] * F[1] * F[5] +
          cste * ds[50] * F[2] * F[3] - cste * ds[41] * F[1] * F[2] +
          2 * s[2] * F[1]) /
        2;
    dP[68] = ((cste * ds[32] * F[4] - 2 * F[0] * ds[14]) * F[8] +
              (2 * ds[14] * F[3] - cste * ds[32] * F[1]) * F[6] -
              cste * ds[50] * F[3] * F[4] + cste * F[0] * ds[50] * F[1]) /
             2;
    dP[77] = ((cste * ds[41] * F[6] - 2 * F[0] * ds[23]) * F[7] -
              cste * ds[50] * F[5] * F[6] - cste * s[5] * F[6] +
              2 * ds[23] * F[4] * F[5] - cste * ds[41] * F[2] * F[4] +
              2 * s[2] * F[4] + cste * F[0] * ds[50] * F[2]) /
             2;
    // (%i21) f90(diff(P[1],F[6]));
    dP[6] = -((2 * ds[6] * F[7] - cste * ds[33] * F[5]) * F[8] -
              cste * ds[42] * F[3] * F[7] + cste * ds[42] * F[1] * F[5] +
              cste * ds[33] * F[2] * F[3] - 2 * ds[6] * F[1] * F[2]) /
            2;
    dP[15] =
        ((cste * ds[33] * F[6] - cste * F[0] * ds[51] + cste * s[3]) * F[7] -
         2 * ds[15] * F[5] * F[6] + cste * ds[51] * F[4] * F[5] -
         2 * s[1] * F[5] - cste * ds[33] * F[2] * F[4] +
         2 * F[0] * ds[15] * F[2]) /
        2;
    dP[24] = ((cste * ds[42] * F[4] - cste * F[0] * ds[51]) * F[8] +
              (cste * ds[51] * F[3] - cste * ds[42] * F[1]) * F[6] -
              2 * ds[24] * F[3] * F[4] + cste * s[5] * F[3] +
              2 * F[0] * ds[24] * F[1] - cste * s[4] * F[1]) /
             2;
    dP[33] = ((2 * ds[6] * F[6] - cste * F[0] * ds[42] + 2 * s[0]) * F[7] -
              cste * ds[33] * F[5] * F[6] + cste * ds[42] * F[4] * F[5] -
              cste * s[3] * F[5] - 2 * ds[6] * F[2] * F[4] +
              cste * F[0] * ds[33] * F[2]) /
             2;
    dP[42] = -((cste * ds[33] * F[7] - 2 * ds[15] * F[5]) * F[8] -
               cste * ds[51] * F[3] * F[7] + cste * ds[51] * F[1] * F[5] +
               2 * ds[15] * F[2] * F[3] - cste * ds[33] * F[1] * F[2]) /
             2;
    dP[51] = ((2 * ds[6] * F[4] - cste * F[0] * ds[33]) * F[8] +
              (cste * ds[33] * F[3] - 2 * ds[6] * F[1]) * F[6] -
              cste * ds[42] * F[3] * F[4] + cste * s[3] * F[3] +
              cste * F[0] * ds[42] * F[1] - 2 * s[0] * F[1]) /
             2;
    dP[60] = -((cste * ds[42] * F[7] - cste * ds[51] * F[5]) * F[8] -
               2 * ds[24] * F[3] * F[7] + 2 * ds[24] * F[1] * F[5] +
               cste * ds[51] * F[2] * F[3] - cste * ds[42] * F[1] * F[2]) /
             2;
    dP[69] = ((cste * ds[33] * F[4] - 2 * F[0] * ds[15]) * F[8] +
              (2 * ds[15] * F[3] - cste * ds[33] * F[1]) * F[6] -
              cste * ds[51] * F[3] * F[4] + 2 * s[1] * F[3] +
              cste * F[0] * ds[51] * F[1] - cste * s[3] * F[1]) /
             2;
    dP[78] = ((cste * ds[42] * F[6] - 2 * F[0] * ds[24] + cste * s[4]) * F[7] -
              cste * ds[51] * F[5] * F[6] + 2 * ds[24] * F[4] * F[5] -
              cste * s[5] * F[5] - cste * ds[42] * F[2] * F[4] +
              cste * F[0] * ds[51] * F[2]) /
             2;
    // (%i22) f90(diff(P[1],F[7]));
    dP[7] = -((2 * ds[7] * F[7] - cste * ds[34] * F[5] + 2 * s[0]) * F[8] -
              cste * ds[43] * F[3] * F[7] + cste * ds[43] * F[1] * F[5] +
              cste * ds[34] * F[2] * F[3] - cste * s[4] * F[3] -
              2 * ds[7] * F[1] * F[2]) /
            2;
    dP[16] = ((cste * ds[34] * F[6] - cste * F[0] * ds[52]) * F[7] -
              2 * ds[16] * F[5] * F[6] + cste * s[3] * F[6] +
              cste * ds[52] * F[4] * F[5] - cste * ds[34] * F[2] * F[4] +
              2 * F[0] * ds[16] * F[2] - cste * F[0] * s[5]) /
             2;
    dP[25] = ((cste * ds[43] * F[4] - cste * F[0] * ds[52]) * F[8] +
              (cste * ds[52] * F[3] - cste * ds[43] * F[1]) * F[6] -
              2 * ds[25] * F[3] * F[4] + 2 * F[0] * ds[25] * F[1]) /
             2;
    dP[34] = ((2 * ds[7] * F[6] - cste * F[0] * ds[43]) * F[7] -
              cste * ds[34] * F[5] * F[6] + 2 * s[0] * F[6] +
              cste * ds[43] * F[4] * F[5] - 2 * ds[7] * F[2] * F[4] +
              cste * F[0] * ds[34] * F[2] - cste * F[0] * s[4]) /
             2;
    dP[43] = -((cste * ds[34] * F[7] - 2 * ds[16] * F[5] + cste * s[3]) * F[8] -
               cste * ds[52] * F[3] * F[7] + cste * ds[52] * F[1] * F[5] +
               2 * ds[16] * F[2] * F[3] - cste * s[5] * F[3] -
               cste * ds[34] * F[1] * F[2]) /
             2;
    dP[52] = ((2 * ds[7] * F[4] - cste * F[0] * ds[34]) * F[8] +
              (cste * ds[34] * F[3] - 2 * ds[7] * F[1]) * F[6] -
              cste * ds[43] * F[3] * F[4] + cste * F[0] * ds[43] * F[1]) /
             2;
    dP[61] =
        -((cste * ds[43] * F[7] - cste * ds[52] * F[5] + cste * s[4]) * F[8] -
          2 * ds[25] * F[3] * F[7] + 2 * ds[25] * F[1] * F[5] +
          cste * ds[52] * F[2] * F[3] - 2 * s[2] * F[3] -
          cste * ds[43] * F[1] * F[2]) /
        2;
    dP[70] = ((cste * ds[34] * F[4] - 2 * F[0] * ds[16]) * F[8] +
              (2 * ds[16] * F[3] - cste * ds[34] * F[1]) * F[6] -
              cste * ds[52] * F[3] * F[4] + cste * F[0] * ds[52] * F[1]) /
             2;
    dP[79] = ((cste * ds[43] * F[6] - 2 * F[0] * ds[25]) * F[7] -
              cste * ds[52] * F[5] * F[6] + cste * s[4] * F[6] +
              2 * ds[25] * F[4] * F[5] - cste * ds[43] * F[2] * F[4] +
              cste * F[0] * ds[52] * F[2] - 2 * F[0] * s[2]) /
             2;
    // (%i23) f90(diff(P[1],F[8]));
    dP[8] = -((2 * ds[8] * F[7] - cste * ds[35] * F[5]) * F[8] -
              cste * ds[44] * F[3] * F[7] + 2 * s[0] * F[7] +
              cste * ds[44] * F[1] * F[5] - cste * s[3] * F[5] +
              cste * ds[35] * F[2] * F[3] - 2 * ds[8] * F[1] * F[2]) /
            2;
    dP[17] = ((cste * ds[35] * F[6] - cste * F[0] * ds[53]) * F[7] -
              2 * ds[17] * F[5] * F[6] + cste * ds[53] * F[4] * F[5] -
              cste * ds[35] * F[2] * F[4] + 2 * F[0] * ds[17] * F[2]) /
             2;
    dP[26] = ((cste * ds[44] * F[4] - cste * F[0] * ds[53]) * F[8] +
              (cste * ds[53] * F[3] - cste * ds[44] * F[1]) * F[6] -
              2 * ds[26] * F[3] * F[4] + cste * s[4] * F[4] +
              2 * F[0] * ds[26] * F[1] - cste * F[0] * s[5]) /
             2;
    dP[35] = ((2 * ds[8] * F[6] - cste * F[0] * ds[44]) * F[7] -
              cste * ds[35] * F[5] * F[6] + cste * ds[44] * F[4] * F[5] -
              2 * ds[8] * F[2] * F[4] + cste * F[0] * ds[35] * F[2]) /
             2;
    dP[44] = -((cste * ds[35] * F[7] - 2 * ds[17] * F[5]) * F[8] -
               cste * ds[53] * F[3] * F[7] + cste * s[3] * F[7] +
               cste * ds[53] * F[1] * F[5] - 2 * s[1] * F[5] +
               2 * ds[17] * F[2] * F[3] - cste * ds[35] * F[1] * F[2]) /
             2;
    dP[53] = ((2 * ds[8] * F[4] - cste * F[0] * ds[35]) * F[8] +
              (cste * ds[35] * F[3] - 2 * ds[8] * F[1]) * F[6] -
              cste * ds[44] * F[3] * F[4] + 2 * s[0] * F[4] +
              cste * F[0] * ds[44] * F[1] - cste * F[0] * s[3]) /
             2;
    dP[62] = -((cste * ds[44] * F[7] - cste * ds[53] * F[5]) * F[8] -
               2 * ds[26] * F[3] * F[7] + cste * s[4] * F[7] +
               2 * ds[26] * F[1] * F[5] - cste * s[5] * F[5] +
               cste * ds[53] * F[2] * F[3] - cste * ds[44] * F[1] * F[2]) /
             2;
    dP[71] = ((cste * ds[35] * F[4] - 2 * F[0] * ds[17]) * F[8] +
              (2 * ds[17] * F[3] - cste * ds[35] * F[1]) * F[6] -
              cste * ds[53] * F[3] * F[4] + cste * s[3] * F[4] +
              cste * F[0] * ds[53] * F[1] - 2 * F[0] * s[1]) /
             2;
    dP[80] = ((cste * ds[44] * F[6] - 2 * F[0] * ds[26]) * F[7] -
              cste * ds[53] * F[5] * F[6] + 2 * ds[26] * F[4] * F[5] -
              cste * ds[44] * F[2] * F[4] + cste * F[0] * ds[53] * F[2]) /
             2;
  }  // end of convertFiniteStrainTangentOperator_PK1_3D

  static void convertFiniteStrainTangentOperator_PK1_2D(
      std::span<mgis::real>& dPs,
      const MaterialDataManager& m,
      const mgis::size_type b,
      const mgis::size_type e) {
    auto* const dP = dPs.data();
    const auto* const ds = m.K.data();
    const auto* const F = m.s1.gradients.data();
    const auto* const s = m.s1.thermodynamic_forces.data();
    // symmetric tensor size
    const auto ss = getStensorSize(Hypothesis::PLANESTRAIN);
    // non symmetric tensor size
    const auto ts = getTensorSize(Hypothesis::PLANESTRAIN);
    // stride associated with Kd
    const auto dP_stride = ts * ts;
    // stride associated with m.K
    const auto ds_stride = ss * ts;
    for (auto i = b; i != e; ++i) {
      auto* const dP_l = dP + dP_stride * i;
      const auto* const ds_l = ds + ds_stride * i;
      const auto* const F_l = F + ts * i;
      const auto* const s_l = s + ss * i;
      convertFiniteStrainTangentOperator_PK1_2D(dP_l, ds_l, F_l, s_l);
    }
  }  // end of convertFiniteStrainTangentOperator_PK1_2D

  static void convertFiniteStrainTangentOperator_PK1_3D(
      std::span<mgis::real>& dPs,
      const MaterialDataManager& m,
      const mgis::size_type b,
      const mgis::size_type e) {
    auto* const dP = dPs.data();
    const auto* const ds = m.K.data();
    const auto* const F = m.s1.gradients.data();
    const auto* const s = m.s1.thermodynamic_forces.data();
    // symmetric tensor size
    const auto ss = getStensorSize(Hypothesis::TRIDIMENSIONAL);
    // non symmetric tensor size
    const auto ts = getTensorSize(Hypothesis::TRIDIMENSIONAL);
    // stride associated with Kd
    const auto dP_stride = ts * ts;
    // stride associated with m.K
    const auto ds_stride = ss * ts;
    for (auto i = b; i != e; ++i) {
      auto* const dP_l = dP + dP_stride * i;
      const auto* const ds_l = ds + ds_stride * i;
      const auto* const F_l = F + ts * i;
      const auto* const s_l = s + ss * i;
      convertFiniteStrainTangentOperator_PK1_3D(dP_l, ds_l, F_l, s_l);
    }
  }  // end of convertFiniteStrainTangentOperator_PK1_3D

  static void convertFiniteStrainTangentOperator_PK1_2D(
      std::span<mgis::real>& K, const MaterialDataManager& m) {
    // check behaviour type
    if (m.b.btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      mgis::raise(
          "convertFiniteStrainTangentOperator: "
          "unsupported behaviour type");
    }
    // non symmetric tensor size
    const auto ts = getTensorSize(m.b.hypothesis);
    // stride associated with K
    const auto Ks = ts * ts;
    // check K size
    if (K.size() != m.n * Ks) {
      mgis::raise(
          "convertFiniteStrainTangentOperator: "
          "unsupported tangent operator");
    }
    convertFiniteStrainTangentOperator_PK1_2D(K, m, 0, m.n);
  }  // end of convertFiniteStrainTangentOperator_PK1_2D

  static void convertFiniteStrainTangentOperator_PK1_3D(
      std::span<mgis::real>& K, const MaterialDataManager& m) {
    // check behaviour type
    if (m.b.btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      mgis::raise(
          "convertFiniteStrainTangentOperator: "
          "unsupported behaviour type");
    }
    // non symmetric tensor size
    const auto ts = getTensorSize(m.b.hypothesis);
    // stride associated with K
    const auto Ks = ts * ts;
    // check K size
    if (K.size() != m.n * Ks) {
      mgis::raise(
          "convertFiniteStrainTangentOperator: "
          "unsupported tangent operator");
    }
    convertFiniteStrainTangentOperator_PK1_3D(K, m, 0, m.n);
  }  // end of convertFiniteStrainTangentOperator_PK1_3D

  void convertFiniteStrainTangentOperator(std::span<mgis::real>& K,
                                          const MaterialDataManager& m,
                                          const FiniteStrainTangentOperator t) {
    const auto h = m.b.hypothesis;
    if (t == FiniteStrainTangentOperator::DPK1_DF) {
      if (h == Hypothesis::TRIDIMENSIONAL) {
        convertFiniteStrainTangentOperator_PK1_3D(K, m);
      } else if ((h == Hypothesis::AXISYMMETRICAL) ||
                 (h == Hypothesis::PLANESTRAIN) ||
                 (h == Hypothesis::GENERALISEDPLANESTRAIN)) {
        convertFiniteStrainTangentOperator_PK1_2D(K, m);
      } else {
        mgis::raise(
            "convertFiniteStrainTangentOperator: "
            "unsupported hypothesis");
      }
    } else {
      mgis::raise(
          "convertFiniteStrainTangentOperator: "
          "unsupported tangent operator");
    }
  }  // end of convertFiniteStrainTangentOperator

  static void convertFiniteStrainTangentOperator_PK1_2D(
      std::span<mgis::real>& K, const BehaviourData& d) {
    // check behaviour type
    if (d.s1.b.btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      mgis::raise(
          "convertFiniteStrainTangentOperator: "
          "unsupported behaviour type");
    }
    // non symmetric tensor size
    const auto ts = getTensorSize(d.s1.b.hypothesis);
    // stride associated with K
    const auto Ks = ts * ts;
    // check K size
    if (K.size() != Ks) {
      mgis::raise(
          "convertFiniteStrainTangentOperator: "
          "unsupported tangent operator");
    }
    convertFiniteStrainTangentOperator_PK1_2D(K.data(), d.K.data(),
                                              d.s1.gradients.data(),
                                              d.s1.thermodynamic_forces.data());
  }  // end of convertFiniteStrainTangentOperator_PK1_2D

  static void convertFiniteStrainTangentOperator_PK1_3D(
      std::span<mgis::real>& K, const BehaviourData& d) {
    // check behaviour type
    if (d.s1.b.btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      mgis::raise(
          "convertFiniteStrainTangentOperator: "
          "unsupported behaviour type");
    }
    // non symmetric tensor size
    const auto ts = getTensorSize(d.s1.b.hypothesis);
    // stride associated with K
    const auto Ks = ts * ts;
    // check K size
    if (K.size() != Ks) {
      mgis::raise(
          "convertFiniteStrainTangentOperator: "
          "unsupported tangent operator");
    }
    convertFiniteStrainTangentOperator_PK1_3D(K.data(), d.K.data(),
                                              d.s1.gradients.data(),
                                              d.s1.thermodynamic_forces.data());
  }  // end of convertFiniteStrainTangentOperator_PK1_3D

  void convertFiniteStrainTangentOperator(std::span<mgis::real>& K,
                                          const BehaviourData& d,
                                          const FiniteStrainTangentOperator t) {
    const auto h = d.s1.b.hypothesis;
    if (t == FiniteStrainTangentOperator::DPK1_DF) {
      if (h == Hypothesis::TRIDIMENSIONAL) {
        convertFiniteStrainTangentOperator_PK1_3D(K, d);
      } else if ((h == Hypothesis::AXISYMMETRICAL) ||
                 (h == Hypothesis::PLANESTRAIN) ||
                 (h == Hypothesis::GENERALISEDPLANESTRAIN)) {
        convertFiniteStrainTangentOperator_PK1_2D(K, d);
      } else {
        mgis::raise(
            "convertFiniteStrainTangentOperator: "
            "unsupported hypothesis");
      }
    } else {
      mgis::raise(
          "convertFiniteStrainTangentOperator: "
          "unsupported tangent operator");
    }
  }  // end of convertFiniteStrainTangentOperator

}  // end of namespace mgis::behaviour
