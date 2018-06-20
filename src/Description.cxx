/*!
 * \file   Description.cxx
 * \brief
 * \author HELFER Thomas 202608
 * \date   20 juin 2018
 */

#include "MFront/Behaviour/Description.hxx"
#include "MFront/LibrariesManager.hxx"
#include "MFront/Raise.hxx"

namespace mfront {

namespace behaviour {

Description::Description() = default;
Description::Description(Description &&) = default;
Description::Description(const Description &) = default;
Description &Description::operator=(Description &&) = default;
Description &Description::operator=(const Description &) = default;
Description::~Description() = default;

Description load(const std::string &l, const std::string &b,
                 const Hypothesis h) {
  auto &lm = mfront::LibrariesManager::get();
  const auto fct = b + '_' + toString(h);
  auto d = Description{};

  d.library = l;
  d.behaviour = b;
  d.function = fct;
  d.hypothesis = h;

  d.tfel_version = lm.getTFELVersion(l, b);
  d.source = lm.getSource(l, b);
  d.btype = [&l, &b, &lm] {
    /* - 0 : general behaviour
     * - 1 : strain based behaviour *
     * - 2 : standard finite strain behaviour *
     * - 3 : cohesive zone model */
    switch (lm.getBehaviourType(l, b)) {
    case 0:
      return Description::GENERALBEHAVIOUR;
      break;
    case 1:
      return Description::STANDARDSTRAINBASEDBEHAVIOUR;
      break;
    case 2:
      return Description::STANDARDFINITESTRAINBEHAVIOUR;
      break;
    case 3:
      return Description::COHESIVEZONEMODEL;
      break;
    }
    raise("load: unsupported behaviour type");
  }();
  // behaviour kinematic
  d.kinematic = [&l, &b, &lm, h] {
    /* - 0: undefined kinematic
     * - 1: standard small strain behaviour kinematic
     * - 2: cohesive zone model kinematic
     * - 3: standard finite strain kinematic (F-Cauchy)
     * - 4: ptest finite strain kinematic (eto-pk1)
     * - 5: Green-Lagrange strain
     * - 6: Miehe Apel Lambrecht logarithmic strain framework */
    switch (lm.getBehaviourKinematic(l, b)) {
    case 0:
      return Description::UNDEFINEDKINEMATIC;
      break;
    case 1:
      return Description::SMALLSTRAINKINEMATIC;
      break;
    case 2:
      return Description::COHESIVEZONEKINEMATIC;
      break;
    case 3:
      return Description::FINITESTRAINKINEMATIC_F_CAUCHY;
      break;
    case 4:
      raise_if(((h != Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN) &&
                (h != Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS)),
               "load: invalid hypothesis for behaviour based on "
               "the eto-pk1 kinematic");
      return Description::FINITESTRAINKINEMATIC_ETO_PK1;
      break;
    }
    raise("load: unsupported behaviour kinematic");
  }();
  // setting the gradients and the fluxes
  switch (d.btype) {
  case Description::GENERALBEHAVIOUR:
    raise("load: general behaviour is not handled yet");
    break;
  case Description::STANDARDSTRAINBASEDBEHAVIOUR:
    raise_if(d.kinematic != Description::SMALLSTRAINKINEMATIC,
             "load: strain based behaviour must be associated with the "
             "small strain kinematic hypothesis");
    d.gradients.push_back({"Strain", Variable::STENSOR});
    d.fluxes.push_back({"Stress", Variable::STENSOR});
    break;
  case Description::COHESIVEZONEMODEL:
    raise_if(d.kinematic != Description::COHESIVEZONEKINEMATIC,
             "load: invalid kinematic assumption for cohesive zone model");
    d.gradients.push_back({"OpeningDisplacement", Variable::VECTOR});
    d.fluxes.push_back({"CohesiveForce", Variable::VECTOR});
    break;
  case Description::STANDARDFINITESTRAINBEHAVIOUR:
    if (d.kinematic == Description::FINITESTRAINKINEMATIC_F_CAUCHY) {
      d.gradients.push_back({"DeformationGradient", Variable::TENSOR});
      d.fluxes.push_back({"CauchyStress", Variable::STENSOR});
    } else if (d.kinematic == Description::FINITESTRAINKINEMATIC_ETO_PK1) {
      raise_if(((h != Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN) &&
                (h != Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS)),
               "load: invalid hypothesis for behaviour based on "
               "the eto-pk1 kinematic");
      d.gradients.push_back({"Strain", Variable::STENSOR});
      d.fluxes.push_back({"PK1-Stress", Variable::TENSOR});
    } else {
      raise("load: invalid kinematic hypothesis for finite strain behaviour");
    }
    break;
  default:
    raise("load: unsupported behaviour type");
  };
  // behaviour symmetry
  d.symmetry = lm.getBehaviourSymmetry(l, b) == 0 ? Description::ISOTROPIC
                                                  : Description::ORTHOTROPIC;
  auto add_mp = [&d](const std::string &mp) {
    d.mps.push_back({mp, Variable::SCALAR});
  };
  if (lm.requiresStiffnessTensor(l, b, h)) {
    if (lm.getElasticStiffnessSymmetry(l, b) == 0) {
      add_mp("YoungModulus");
      add_mp("PoissonRatio");
    } else {
      raise_if(d.symmetry != Description::ORTHOTROPIC,
               "load: the behaviour must be orthotropic "
               "for the elastic stiffness symmetry to be orthotropic");
      add_mp("YoungModulus1");
      add_mp("YoungModulus2");
      add_mp("YoungModulus3");
      add_mp("PoissonRatio12");
      add_mp("PoissonRatio23");
      add_mp("PoissonRatio13");
      if ((h == Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN) ||
          (h == Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS)) {
      } else if ((h == Hypothesis::PLANESTRESS) ||
                 (h == Hypothesis::PLANESTRAIN) ||
                 (h == Hypothesis::AXISYMMETRICAL) ||
                 (h == Hypothesis::GENERALISEDPLANESTRAIN)) {
        add_mp("ShearModulus12");
      } else if (h == Hypothesis::TRIDIMENSIONAL) {
        add_mp("ShearModulus12");
        add_mp("ShearModulus23");
        add_mp("ShearModulus13");
      }
    }
  }
  if (lm.requiresThermalExpansionCoefficientTensor(l, b, h)) {
    if (d.symmetry == Description::ORTHOTROPIC) {
      add_mp("ThermalExpansion1");
      add_mp("ThermalExpansion2");
      add_mp("ThermalExpansion3");
    } else {
      add_mp("ThermalExpansion");
    }
  }
  // standard material properties
  for (const auto &mp : lm.getMaterialPropertiesNames(l, b, h)) {
    add_mp(mp);
  }
  // internal state variables
  const auto ivnames = lm.getInternalStateVariablesNames(l, b, h);
  const auto ivtypes = lm.getInternalStateVariablesTypes(l, b, h);
  raise_if(ivnames.size() != ivtypes.size(),
           "load: the number of internal state variables names does not match "
           "the number of internal state variables types");
  for (decltype(ivnames.size()) i = 0; i != ivnames.size(); ++i) {
    switch (ivtypes[i]) {
    case 0:
      d.isvs.push_back({ivnames[i], Variable::SCALAR});
      break;
    default:
      raise("load: unsupported internal state variable type");
    }
  }
  // external state variables
  for (const auto &esv : lm.getExternalStateVariablesNames(l, b, h)) {
    d.esvs.push_back({esv, Variable::SCALAR});
  }
  //       //! parameters
  //       const auto pn = lm.getUMATParametersNames(l, fct, h);
  //       const auto pt = lm.getUMATParametersTypes(l, fct, h);
  //       throw_if(
  //           pn.size() != pt.size(),
  //           "inconsistent size between parameters' names and
  //           parameters' sizes");
  //       for (decltype(pn.size()) i = 0; i != pn.size(); ++i) {
  //         if (pt[i] == 0) {
  //           d.pnames.push_back(pn[i]);
  //         } else if (pt[i] == 1) {
  //           d.ipnames.push_back(pn[i]);
  //         } else if (pt[i] == 2) {
  //           d.upnames.push_back(pn[i]);
  //         } else {
  //           throw_if(true,
  //                    "unsupported parameter type for parameter '" +
  //                    pn[i]
  //                    + "'");
  //         }
  //       }
  return d;
} // end of load

} // end of namespace behaviour

} // end of namespace mfront
