/*!
 * \file   src/Behaviour.cxx
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

#include <cstdlib>
#include <iterator>
#include "MGIS/Raise.hxx"
#include "MGIS/Context.hxx"
#include "MGIS/LibrariesManager.hxx"
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"

namespace mgis::behaviour {

  static std::array<mgis::real, 9u> buildRotationMatrix(const real *const a) {
    return buildRotationMatrix(std::span<const mgis::real, 2u>{a, 2u});
  }  // end of buildRotationMatrix

  static std::array<mgis::real, 9u> buildRotationMatrix(const real *const a1,
                                                        const real *const a2) {
    return buildRotationMatrix(std::span<const mgis::real, 3u>{a1, 3u},
                               std::span<const mgis::real, 3u>{a2, 3u});
  }  // end of buildRotationMatrix

  BehaviourInitializeFunction::BehaviourInitializeFunction() = default;
  BehaviourInitializeFunction::BehaviourInitializeFunction(
      BehaviourInitializeFunction &&) = default;
  BehaviourInitializeFunction::BehaviourInitializeFunction(
      const BehaviourInitializeFunction &) = default;
  BehaviourInitializeFunction &BehaviourInitializeFunction::operator=(
      BehaviourInitializeFunction &&) = default;
  BehaviourInitializeFunction &BehaviourInitializeFunction::operator=(
      const BehaviourInitializeFunction &) = default;
  BehaviourInitializeFunction::~BehaviourInitializeFunction() = default;

  BehaviourPostProcessing::BehaviourPostProcessing() = default;
  BehaviourPostProcessing::BehaviourPostProcessing(BehaviourPostProcessing &&) =
      default;
  BehaviourPostProcessing::BehaviourPostProcessing(
      const BehaviourPostProcessing &) = default;
  BehaviourPostProcessing &BehaviourPostProcessing::operator=(
      BehaviourPostProcessing &&) = default;
  BehaviourPostProcessing &BehaviourPostProcessing::operator=(
      const BehaviourPostProcessing &) = default;
  BehaviourPostProcessing::~BehaviourPostProcessing() = default;

  Behaviour::Behaviour() = default;
  Behaviour::Behaviour(Behaviour &&) = default;
  Behaviour::Behaviour(const Behaviour &) = default;
  Behaviour &Behaviour::operator=(Behaviour &&) = default;
  Behaviour &Behaviour::operator=(const Behaviour &) = default;
  Behaviour::~Behaviour() = default;

  static Behaviour load_behaviour(const std::string &l,
                                  const std::string &b,
                                  const Hypothesis h) {
    auto &lm = mgis::LibrariesManager::get();
    //
    auto d = Behaviour{};
    loadBehaviourDescription(d, l, b, h);
    d.b = lm.getBehaviour(l, b, h);
    //
    if (d.btype == BehaviourDescription::STANDARDFINITESTRAINBEHAVIOUR) {
      d.options.resize(2, mgis::real(0));
    }
    // initialize functions
    for (const auto &i : lm.getBehaviourInitializeFunctions(l, b, h)) {
      BehaviourInitializeFunction ifct;
      ifct.f = lm.getBehaviourInitializeFunction(l, b, i, h);
      ifct.inputs = getBehaviourInitializeFunctionInputs(l, b, i, h);
      d.initialize_functions.insert({i, ifct});
    }
    // post-processings
    for (const auto &i : lm.getBehaviourPostProcessings(l, b, h)) {
      BehaviourPostProcessing pfct;
      pfct.f = lm.getBehaviourPostProcessing(l, b, i, h);
      pfct.outputs = getBehaviourPostProcessingOutputs(l, b, i, h);
      d.postprocessings.insert({i, pfct});
    }
    return d;
  }  // end of load_behaviour

  Behaviour load(const std::string &l,
                 const std::string &b,
                 const Hypothesis h) {
    if (isStandardFiniteStrainBehaviour(l, b)) {
      mgis::raise(
          "mgis::behaviour::load: "
          "This version of the load function shall not be called "
          "for finite strain behaviour: you shall specify finite "
          "strain options");
    }
    auto d = load_behaviour(l, b, h);
    if (d.symmetry == Behaviour::ORTHOTROPIC) {
      auto &lm = mgis::LibrariesManager::get();
      d.rotate_gradients_ptr = lm.getRotateBehaviourGradientsFunction(l, b, h);
      d.rotate_array_of_gradients_ptr =
          lm.getRotateArrayOfBehaviourGradientsFunction(l, b, h);
      d.rotate_thermodynamic_forces_ptr =
          lm.getRotateBehaviourThermodynamicForcesFunction(l, b, h);
      d.rotate_array_of_thermodynamic_forces_ptr =
          lm.getRotateArrayOfBehaviourThermodynamicForcesFunction(l, b, h);
      d.rotate_tangent_operator_blocks_ptr =
          lm.getRotateBehaviourTangentOperatorBlocksFunction(l, b, h);
      d.rotate_array_of_tangent_operator_blocks_ptr =
          lm.getRotateArrayOfBehaviourTangentOperatorBlocksFunction(l, b, h);
    }
    return d;
  }  // end of load

  Behaviour load(const FiniteStrainBehaviourOptions &o,
                 const std::string &l,
                 const std::string &b,
                 const Hypothesis h) {
    auto d = load_behaviour(l, b, h);
    if (d.btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      mgis::raise(
          "mgis::behaviour::load: "
          "This method shall only be called for finite strain behaviour");
    }
    if (o.stress_measure == FiniteStrainBehaviourOptions::CAUCHY) {
      d.options[0] = mgis::real(0);
    } else if (o.stress_measure == FiniteStrainBehaviourOptions::PK2) {
      d.options[0] = mgis::real(1);
      d.thermodynamic_forces[0] = {"SecondPiolaKirchhoffStress",
                                   Variable::STENSOR, 1};
    } else if (o.stress_measure == FiniteStrainBehaviourOptions::PK1) {
      d.options[0] = mgis::real(2);
      d.thermodynamic_forces[0] = {"FirstPiolaKirchhoffStress",
                                   Variable::TENSOR, 3};
    } else {
      mgis::raise(
          "mgis::behaviour::load: "
          "internal error (unsupported stress measure)");
    }
    if (o.tangent_operator == FiniteStrainBehaviourOptions::DSIG_DF) {
      d.options[1] = mgis::real(0);
    } else if (o.tangent_operator == FiniteStrainBehaviourOptions::DS_DEGL) {
      d.options[1] = mgis::real(1);
      d.to_blocks[0] = {{"SecondPiolaKirchhoffStress", Variable::STENSOR, 1},
                        {"GreenLagrangeStrain", Variable::STENSOR, 1}};

    } else if (o.tangent_operator == FiniteStrainBehaviourOptions::DPK1_DF) {
      d.options[1] = mgis::real(2);
      d.to_blocks[0] = {{"FirstPiolaKirchhoffStress", Variable::TENSOR, 3},
                        {"DeformationGradient", Variable::TENSOR, 3}};
    } else if (o.tangent_operator == FiniteStrainBehaviourOptions::DTAU_DDF) {
      d.options[1] = mgis::real(3);
      d.to_blocks[0] = {
          {"KirchhoffStress", Variable::STENSOR, 3},
          {"SpatialIncrementOfTheDeformationGradient", Variable::TENSOR, 3}};
    } else {
      mgis::raise(
          "mgis::behaviour::load: "
          "internal error (unsupported tangent operator)");
    }
    if (d.symmetry == Behaviour::ORTHOTROPIC) {
      auto &lm = mgis::LibrariesManager::get();
      d.rotate_gradients_ptr = lm.getRotateBehaviourGradientsFunction(l, b, h);
      d.rotate_array_of_gradients_ptr =
          lm.getRotateArrayOfBehaviourGradientsFunction(l, b, h);
      d.rotate_thermodynamic_forces_ptr =
          lm.getRotateBehaviourThermodynamicForcesFunction(l, b, h,
                                                           o.stress_measure);
      d.rotate_array_of_thermodynamic_forces_ptr =
          lm.getRotateArrayOfBehaviourThermodynamicForcesFunction(
              l, b, h, o.stress_measure);
      d.rotate_tangent_operator_blocks_ptr =
          lm.getRotateBehaviourTangentOperatorBlocksFunction(
              l, b, h, o.tangent_operator);
      d.rotate_array_of_tangent_operator_blocks_ptr =
          lm.getRotateArrayOfBehaviourTangentOperatorBlocksFunction(
              l, b, h, o.tangent_operator);
    }
    return d;
  }  // end of load

  std::optional<Behaviour> load(Context &ctx,
                                const std::string &l,
                                const std::string &f,
                                const Hypothesis h) noexcept {
    try {
      return load(l, f, h);
    } catch (...) {
      registerExceptionInErrorBacktrace(ctx);
    }
    return {};
  } // end of load

  std::optional<Behaviour> load(Context &ctx,
                                const FiniteStrainBehaviourOptions &o,
                                const std::string &l,
                                const std::string &f,
                                const Hypothesis h) noexcept {
    try {
      return load(o, l, f, h);
    } catch (...) {
      registerExceptionInErrorBacktrace(ctx);
    }
    return {};
  }

  void rotateGradients(std::span<real> g,
                       const Behaviour &b,
                       const std::span<const real> &r) {
    rotateGradients(g, b, g, r);
  }  // end of rotateGradients

  void rotateGradients(std::span<real> g,
                       const Behaviour &b,
                       const RotationMatrix2D &r) {
    rotateGradients(g, b, g, r);
  }  // end of rotateGradients

  void rotateGradients(std::span<real> g,
                       const Behaviour &b,
                       const RotationMatrix3D &r) {
    rotateGradients(g, b, g, r);
  }  // end of rotateGradients

  static mgis::size_type checkRotateFunctionInputs(
      const char *const m,
      const std::span<const real> &mv,
      const std::span<const real> &gv,
      const mgis::size_type vsize) {
    if (gv.size() == 0) {
      mgis::raise(std::string(m) + ": no values given for the gradients");
    }
    auto dv = std::div(static_cast<long long>(gv.size()), vsize);
    if (dv.rem != 0) {
      mgis::raise(std::string(m) +
                  ": invalid array size in the global frame "
                  "(not a multiple of the gradients size)");
    }
    if (mv.size() != gv.size()) {
      mgis::raise(std::string(m) + ": unmatched array sizes");
    }
    return dv.quot;
  }

  static void checkRotationMatrix2D(const char *const m,
                                    const RotationMatrix2D &r,
                                    const Behaviour &b,
                                    const mgis::size_type nipts) {
    if (getSpaceDimension(b.hypothesis) != 2u) {
      mgis::raise(std::string(m) + ": a 2D rotation matrix can't be used in '" +
                  toString(b.hypothesis) + "'");
    }
    if (r.a.size() != 2u) {
      if (nipts != r.a.size() % 2) {
        mgis::raise(std::string(m) +
                    ": the number of integration points handled "
                    "by the rotation matrix is different from the "
                    "number of integration points of the field to be rotated");
      }
    }
  }  // end of checkRotationMatrix2D

  static void checkRotationMatrix3D(const char *const m,
                                    const RotationMatrix3D &r,
                                    const Behaviour &b,
                                    const mgis::size_type nipts) {
    if (getSpaceDimension(b.hypothesis) != 3u) {
      mgis::raise(std::string(m) + ": a 3D rotation matrix can't be used in '" +
                  toString(b.hypothesis) + "'");
    }
    if (r.a1.a.size() != 3u) {
      if (nipts != r.a1.a.size() % 3) {
        mgis::raise(std::string(m) +
                    ": the number of integration points handled "
                    "by the rotation matrix is different from the "
                    "number of integration points of the field to be rotated");
      }
    }
    if (r.a2.a.size() != 3u) {
      if (nipts != r.a2.a.size() % 3) {
        mgis::raise(std::string(m) +
                    ": the number of integration points handled "
                    "by the rotation matrix is different from the "
                    "number of integration points of the field to be rotated");
      }
    }
  }  // end of checkRotationMatrix3D

  static void checkBehaviourRotateGradients(const Behaviour &b) {
    if ((b.rotate_gradients_ptr == nullptr) ||
        (b.rotate_array_of_gradients_ptr == nullptr)) {
      mgis::raise(
          "rotateGradients: no function performing the rotation of "
          "the gradients defined");
    }
  }  // end of checkBehaviourRotateGradients

  void rotateGradients(std::span<real> mg,
                       const Behaviour &b,
                       const std::span<const real> &gg,
                       const std::span<const real> &r) {
    checkBehaviourRotateGradients(b);
    const auto gsize = getArraySize(b.gradients, b.hypothesis);
    const auto nipts =
        checkRotateFunctionInputs("rotateGradients", mg, gg, gsize);
    if (r.size() == 0) {
      mgis::raise("rotateGradients: no values given for the rotation matrices");
    }
    const auto rdv = std::div(static_cast<long long>(r.size()), size_type{9});
    if (rdv.rem != 0) {
      mgis::raise(
          "rotateGradients: "
          "invalid size for the rotation matrix array");
    }
    if (rdv.quot == 1) {
      b.rotate_array_of_gradients_ptr(mg.data(), gg.data(), r.data(), nipts);
    } else {
      if (rdv.quot != nipts) {
        mgis::raise(
            "the number of integration points for the gradients does not match "
            "the number of integration points for the rotation matrices (" +
            std::to_string(nipts) + " vs " + std::to_string(rdv.quot) + ")");
      }
      for (size_type i = 0; i != nipts; ++i) {
        const auto o = i * gsize;
        b.rotate_gradients_ptr(mg.data() + o, gg.data() + o, r.data() + i * 9);
      }
    }
  }  // end of rotateGradients

  void rotateGradients(std::span<real> mg,
                       const Behaviour &b,
                       const std::span<const real> &gg,
                       const RotationMatrix2D &r) {
    checkBehaviourRotateGradients(b);
    const auto gsize = getArraySize(b.gradients, b.hypothesis);
    const auto nipts =
        checkRotateFunctionInputs("rotateGradients", mg, gg, gsize);
    checkRotationMatrix2D("rotateGradients", r, b, nipts);
    if (r.a.size() == 2u) {
      const auto m = buildRotationMatrix(r.a.data());
      b.rotate_array_of_gradients_ptr(mg.data(), gg.data(), m.data(), nipts);
    } else {
      for (size_type i = 0; i != nipts; ++i) {
        const auto m = buildRotationMatrix(r.a.data() + 2 * i);
        const auto o = i * gsize;
        b.rotate_gradients_ptr(mg.data() + o, gg.data() + o, m.data());
      }
    }
  }  // end of rotateGradients

  void rotateGradients(std::span<real> mg,
                       const Behaviour &b,
                       const std::span<const real> &gg,
                       const RotationMatrix3D &r) {
    checkBehaviourRotateGradients(b);
    const auto gsize = getArraySize(b.gradients, b.hypothesis);
    const auto nipts =
        checkRotateFunctionInputs("rotateGradients", mg, gg, gsize);
    checkRotationMatrix3D("rotateGradients", r, b, nipts);
    if ((r.a1.a.size() == 3u) && ((r.a2.a.size() == 3u))) {
      const auto m = buildRotationMatrix(r.a1.a.data(), r.a2.a.data());
      b.rotate_array_of_gradients_ptr(mg.data(), gg.data(), m.data(), nipts);
    } else {
      const auto o1 = (r.a1.a.size() == 3u) ? 0u : 3u;
      const auto o2 = (r.a2.a.size() == 3u) ? 0u : 3u;
      for (size_type i = 0; i != nipts; ++i) {
        const auto m = buildRotationMatrix(r.a1.a.data() + o1 * i,  //
                                           r.a2.a.data() + o2 * i);
        const auto o = i * gsize;
        b.rotate_gradients_ptr(mg.data() + o, gg.data() + o, m.data());
      }
    }
  }  // end of rotateGradients

  void rotateThermodynamicForces(std::span<real> tf,
                                 const Behaviour &b,
                                 const std::span<const real> &r) {
    rotateThermodynamicForces(tf, b, tf, r);
  }  // end of rotateThermodynamicForces

  void rotateThermodynamicForces(std::span<real> tf,
                                 const Behaviour &b,
                                 const RotationMatrix2D &r) {
    rotateThermodynamicForces(tf, b, tf, r);
  }  // end of rotateThermodynamicForces

  void rotateThermodynamicForces(std::span<real> tf,
                                 const Behaviour &b,
                                 const RotationMatrix3D &r) {
    rotateThermodynamicForces(tf, b, tf, r);
  }  // end of rotateThermodynamicForces

  static void checkBehaviourRotateThermodynamicForces(const Behaviour &b) {
    if ((b.rotate_thermodynamic_forces_ptr == nullptr) ||
        (b.rotate_array_of_thermodynamic_forces_ptr == nullptr)) {
      mgis::raise(
          "rotateThermodynamicForces: no function performing the rotation of "
          "the thermodynamic forces defined");
    }
  }  // end of checkBehaviourRotateThermodynamicForces

  void rotateThermodynamicForces(std::span<real> gtf,
                                 const Behaviour &b,
                                 const std::span<const real> &mtf,
                                 const std::span<const real> &r) {
    checkBehaviourRotateThermodynamicForces(b);
    const auto tfsize = getArraySize(b.thermodynamic_forces, b.hypothesis);
    const auto nipts = checkRotateFunctionInputs("rotateThermodynamicForces",
                                                 mtf, gtf, tfsize);
    if (r.size() == 0) {
      mgis::raise(
          "rotateThermodynamicForces: "
          "no values given for the rotation matrices");
    }
    auto rdv = std::div(static_cast<long long>(r.size()), size_type{9});
    if (rdv.rem != 0) {
      mgis::raise(
          "rotateThermodynamicForces: "
          "invalid size for the rotation matrix array");
    }
    if (rdv.quot == 1) {
      b.rotate_array_of_thermodynamic_forces_ptr(gtf.data(), mtf.data(),
                                                 r.data(), nipts);
    } else {
      if (rdv.quot != nipts) {
        mgis::raise(
            "rotateThermodynamicForces: "
            "the number of integration points for the thermodynamic forces "
            "does not match the number of integration points for the rotation "
            "matrices (" +
            std::to_string(nipts) + " vs " + std::to_string(rdv.quot) + ")");
      }
      for (size_type i = 0; i != nipts; ++i) {
        const auto o = i * tfsize;
        b.rotate_thermodynamic_forces_ptr(gtf.data() + o, mtf.data() + o,
                                          r.data() + i * 9);
      }
    }
  }  // end of rotateThermodynamicForces

  void rotateThermodynamicForces(std::span<real> gtf,
                                 const Behaviour &b,
                                 const std::span<const real> &mtf,
                                 const RotationMatrix2D &r) {
    checkBehaviourRotateThermodynamicForces(b);
    const auto tfsize = getArraySize(b.thermodynamic_forces, b.hypothesis);
    const auto nipts = checkRotateFunctionInputs("rotateThermodynamicForces",
                                                 gtf, mtf, tfsize);
    checkRotationMatrix2D("rotateThermodynamicForces", r, b, nipts);
    if (r.a.size() == 2u) {
      const auto m = buildRotationMatrix(r.a.data());
      b.rotate_array_of_thermodynamic_forces_ptr(gtf.data(), mtf.data(),
                                                 m.data(), nipts);
    } else {
      for (size_type i = 0; i != nipts; ++i) {
        const auto m = buildRotationMatrix(r.a.data() + 2 * i);
        const auto o = i * tfsize;
        b.rotate_thermodynamic_forces_ptr(gtf.data() + o, mtf.data() + o,
                                          m.data());
      }
    }
  }  // end of rotateThermodynamicForces

  void rotateThermodynamicForces(std::span<real> gtf,
                                 const Behaviour &b,
                                 const std::span<const real> &mtf,
                                 const RotationMatrix3D &r) {
    checkBehaviourRotateThermodynamicForces(b);
    const auto tfsize = getArraySize(b.thermodynamic_forces, b.hypothesis);
    const auto nipts = checkRotateFunctionInputs("rotateThermodynamicForces",
                                                 gtf, mtf, tfsize);
    checkRotationMatrix3D("rotateThermodynamicForces", r, b, nipts);
    if ((r.a1.a.size() == 3u) && ((r.a2.a.size() == 3u))) {
      const auto m = buildRotationMatrix(r.a1.a.data(), r.a2.a.data());
      b.rotate_array_of_thermodynamic_forces_ptr(gtf.data(), mtf.data(),
                                                 m.data(), nipts);
    } else {
      const auto o1 = (r.a1.a.size() == 3u) ? 0u : 3u;
      const auto o2 = (r.a2.a.size() == 3u) ? 0u : 3u;
      for (size_type i = 0; i != nipts; ++i) {
        const auto m = buildRotationMatrix(r.a1.a.data() + o1 * i,  //
                                           r.a2.a.data() + o2 * i);
        const auto o = i * tfsize;
        b.rotate_thermodynamic_forces_ptr(gtf.data() + o, mtf.data() + o,
                                          m.data());
      }
    }
  }  // end of rotateThermodynamicForces

  void rotateTangentOperatorBlocks(std::span<real> K,
                                   const Behaviour &b,
                                   const std::span<const real> &r) {
    rotateTangentOperatorBlocks(K, b, K, r);
  }  // end of rotateTangentOperatorBlocks

  void rotateTangentOperatorBlocks(std::span<real> K,
                                   const Behaviour &b,
                                   const RotationMatrix2D &r) {
    rotateTangentOperatorBlocks(K, b, K, r);
  }  // end of rotateTangentOperatorBlocks

  void rotateTangentOperatorBlocks(std::span<real> K,
                                   const Behaviour &b,
                                   const RotationMatrix3D &r) {
    rotateTangentOperatorBlocks(K, b, K, r);
  }  // end of rotateTangentOperatorBlocks

  static void checkBehaviourRotateTangentOperatorBlocks(const Behaviour &b) {
    if ((b.rotate_tangent_operator_blocks_ptr == nullptr) ||
        (b.rotate_array_of_tangent_operator_blocks_ptr == nullptr)) {
      mgis::raise(
          "rotateTangentOperatorBlocks: no function performing the rotation of "
          "the thermodynamic forces defined");
    }
  }  // end of checkBehaviourRotateTangentOperatorBlocks

  void rotateTangentOperatorBlocks(std::span<real> gK,
                                   const Behaviour &b,
                                   const std::span<const real> &mK,
                                   const std::span<const real> &r) {
    checkBehaviourRotateTangentOperatorBlocks(b);
    const auto Ksize = getTangentOperatorArraySize(b);
    const auto nipts =
        checkRotateFunctionInputs("rotateTangentOperatorBlocks", mK, gK, Ksize);
    if (r.size() == 0) {
      mgis::raise(
          "rotateTangentOperatorBlocks: "
          "empty array for the rotation matrix");
    }
    if (mK.size() != gK.size()) {
      mgis::raise("rotateTangentOperatorBlocks: unmatched array sizes");
    }
    auto rdv = std::div(static_cast<long long>(r.size()), size_type{9});
    if (rdv.rem != 0) {
      mgis::raise(
          "rotateTangentOperatorBlocks: "
          "invalid size for the rotation matrix array");
    }
    if (rdv.quot == 1) {
      b.rotate_array_of_tangent_operator_blocks_ptr(gK.data(), mK.data(),
                                                    r.data(), nipts);
    } else {
      if (rdv.quot != nipts) {
        mgis::raise(
            "the number of integration points for the tangent operators does "
            "not match the number of integration points for the rotation "
            "matrices (" +
            std::to_string(nipts) + " vs " + std::to_string(rdv.quot) + ")");
      }
      for (size_type i = 0; i != nipts; ++i) {
        const auto o = i * Ksize;
        b.rotate_tangent_operator_blocks_ptr(gK.data() + o, mK.data() + o,
                                             r.data() + 9 * i);
      }
    }
  }  // end of rotateTangentOperatorBlocks

  void rotateTangentOperatorBlocks(std::span<real> gK,
                                   const Behaviour &b,
                                   const std::span<const real> &mK,
                                   const RotationMatrix2D &r) {
    checkBehaviourRotateTangentOperatorBlocks(b);
    const auto Ksize = getTangentOperatorArraySize(b);
    const auto nipts =
        checkRotateFunctionInputs("rotateTangentOperatorBlocks", gK, mK, Ksize);
    checkRotationMatrix2D("rotateTangentOperatorBlocks", r, b, nipts);
    if (r.a.size() == 2u) {
      const auto m = buildRotationMatrix(r.a.data());
      b.rotate_array_of_tangent_operator_blocks_ptr(gK.data(), mK.data(),
                                                    m.data(), nipts);
    } else {
      for (size_type i = 0; i != nipts; ++i) {
        const auto m = buildRotationMatrix(r.a.data() + 2 * i);
        const auto o = i * Ksize;
        b.rotate_tangent_operator_blocks_ptr(gK.data() + o, mK.data() + o,
                                             m.data());
      }
    }
  }  // end of rotateTangentOperatorBlocks

  void rotateTangentOperatorBlocks(std::span<real> gK,
                                   const Behaviour &b,
                                   const std::span<const real> &mK,
                                   const RotationMatrix3D &r) {
    checkBehaviourRotateTangentOperatorBlocks(b);
    const auto Ksize = getTangentOperatorArraySize(b);
    const auto nipts =
        checkRotateFunctionInputs("rotateTangentOperatorBlocks", gK, mK, Ksize);
    checkRotationMatrix3D("rotateTangentOperatorBlocks", r, b, nipts);
    if ((r.a1.a.size() == 3u) && ((r.a2.a.size() == 3u))) {
      const auto m = buildRotationMatrix(r.a1.a.data(), r.a2.a.data());
      b.rotate_array_of_tangent_operator_blocks_ptr(gK.data(), mK.data(),
                                                    m.data(), nipts);
    } else {
      const auto o1 = (r.a1.a.size() == 3u) ? 0u : 3u;
      const auto o2 = (r.a2.a.size() == 3u) ? 0u : 3u;
      for (size_type i = 0; i != nipts; ++i) {
        const auto m = buildRotationMatrix(r.a1.a.data() + o1 * i,  //
                                           r.a2.a.data() + o2 * i);
        const auto o = i * Ksize;
        b.rotate_tangent_operator_blocks_ptr(gK.data() + o, mK.data() + o,
                                             m.data());
      }
    }
  }  // end of rotateTangentOperatorBlocks

  void setParameter(const Behaviour &b, const std::string &n, const double v) {
    auto &lm = mgis::LibrariesManager::get();
    lm.setParameter(b.library, b.behaviour, b.hypothesis, n, v);
  }  // end of setParameter

  void setParameter(const Behaviour &b, const std::string &n, const int v) {
    auto &lm = mgis::LibrariesManager::get();
    lm.setParameter(b.library, b.behaviour, b.hypothesis, n, v);
  }  // end of setParameter

  void setParameter(const Behaviour &b,
                    const std::string &n,
                    const unsigned short v) {
    auto &lm = mgis::LibrariesManager::get();
    lm.setParameter(b.library, b.behaviour, b.hypothesis, n, v);
  }  // end of setParameter

  size_type getInitializeFunctionVariablesArraySize(const Behaviour &b,
                                                    const std::string_view n) {
    const auto p = b.initialize_functions.find(n);
    if (p == b.initialize_functions.end()) {
      mgis::raise(
          "getInitializeFunctionVariables: "
          "no initialize function named '" +
          std::string{n} + "'");
    }
    return getArraySize(p->second.inputs, b.hypothesis);
  }  // end of getInitializeFunctionVariablesArraySize

  std::vector<mgis::real> allocateInitializeFunctionVariables(
      const Behaviour &b, const std::string_view n) {
    const auto s = getInitializeFunctionVariablesArraySize(b, n);
    std::vector<mgis::real> inputs;
    inputs.resize(s, real{0});
    return inputs;
  }  // end of allocateInitializeFunctionVariables

  size_type getPostProcessingVariablesArraySize(const Behaviour &b,
                                                const std::string_view n) {
    const auto p = b.postprocessings.find(n);
    if (p == b.postprocessings.end()) {
      mgis::raise(
          "getPostProcessingVariables: "
          "no post-processing named '" +
          std::string{n} + "'");
    }
    return getArraySize(p->second.outputs, b.hypothesis);
  }  // end of getPostProcessingVariablesArraySize

  std::vector<mgis::real> allocatePostProcessingVariables(
      const Behaviour &b, const std::string_view n) {
    const auto s = getPostProcessingVariablesArraySize(b, n);
    std::vector<mgis::real> outputs;
    outputs.resize(s, real{0});
    return outputs;
  }

}  // end of namespace mgis::behaviour
