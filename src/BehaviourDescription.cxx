/*!
 * \file   src/BehaviourDescription.cxx
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
#include "MGIS/LibrariesManager.hxx"
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MGIS/Behaviour/BehaviourDescription.hxx"

namespace mgis::behaviour {

  template <typename ErrorHandler>
  static std::vector<Variable> buildVariablesList(
      ErrorHandler &raise,
      const std::vector<std::string> &names,
      const std::vector<int> &types) {
    std::vector<Variable> vars;
    if (names.size() != types.size()) {
      raise(
          "the number of internal state variables names does not match "
          "the number of internal state variables types");
    }
    for (decltype(names.size()) i = 0; i != names.size(); ++i) {
      vars.push_back({names[i], getVariableType(types[i]), types[i]});
    }
    return vars;
  }  // end of buildVariablesList

  template <typename ErrorHandler>
  static void checkGradientsAndThermodynamicForcesConsistency(
      ErrorHandler &raise,
      const std::vector<Variable> &gradients,
      const std::vector<Variable> &thermodynamic_forces,
      const Variable &g,
      const Variable &t) {
    auto raise_if = [&raise](const bool c, const std::string &m) {
      if (c) {
        raise(m);
      }
    };
    raise_if(gradients.size() != thermodynamic_forces.size(),
             "the number of gradients does not match the number of "
             "thermodynamic forces");
    raise_if(gradients.size() != 1u, "invalid number of gradients");
    raise_if(gradients[0].name != g.name, "invalid gradient name");
    raise_if(gradients[0].type != g.type, "invalid gradient type");
    raise_if(thermodynamic_forces[0].name != t.name,
             "invalid thermodynamic force name");
    raise_if(thermodynamic_forces[0].type != t.type,
             "invalid thermodynamic force type");
  }  // end of checkGradientsAndThermodynamicForcesConsistency

  static std::pair<Variable, Variable> getJacobianBlockVariables(
      const BehaviourDescription &b,
      const std::pair<std::string, std::string> &block) {
    auto found = false;
    std::pair<Variable, Variable> v;
    auto assign_if = [&found, &block, &v](const Variable &v1,
                                          const Variable &v2) {
      mgis::raise_if(found,
                     "getJacobianBlockVariables: "
                     "multiple definition for block {" +
                         block.first + "," + block.second + "}");
      found = true;
      v = {v1, v2};
    };
    for (const auto &f : b.thermodynamic_forces) {
      for (const auto &g : b.gradients) {
        if ((block.first == f.name) && (block.second == g.name)) {
          assign_if(f, g);
        }
      }
      for (const auto &e : b.esvs) {
        if ((block.first == f.name) && (block.second == e.name)) {
          assign_if(f, e);
        }
      }
    }
    for (const auto &i : b.isvs) {
      for (const auto &g : b.gradients) {
        if ((block.first == i.name) && (block.second == g.name)) {
          assign_if(i, g);
        }
      }
      for (const auto &e : b.esvs) {
        if ((block.first == i.name) && (block.second == e.name)) {
          assign_if(i, e);
        }
      }
    }
    if (!found) {
      mgis::raise(
          "getJacobianBlockVariables: "
          "tangent operator block {" +
          block.first + "," + block.second + "} is invalid");
    }
    return v;
  }  // end of getJacobianBlockVariables

  BehaviourDescription::BehaviourDescription() = default;
  BehaviourDescription::BehaviourDescription(BehaviourDescription &&) = default;
  BehaviourDescription::BehaviourDescription(const BehaviourDescription &) =
      default;
  BehaviourDescription &BehaviourDescription::operator=(
      BehaviourDescription &&) = default;
  BehaviourDescription &BehaviourDescription::operator=(
      const BehaviourDescription &) = default;
  BehaviourDescription::~BehaviourDescription() = default;

  void loadBehaviourDescription(BehaviourDescription &d,
                                const std::string &l,
                                const std::string &b,
                                const Hypothesis h) {
    auto &lm = mgis::LibrariesManager::get();
    const auto fct = b + '_' + toString(h);
    auto raise = [&b, &l](const std::string &msg) {
      mgis::raise("load: " + msg + ".\nError while trying to load behaviour '" +
                  b + "' in library '" + l + "'\n");
    };
    auto raise_if = [&b, &l](const bool c, const std::string &msg) {
      if (c) {
        mgis::raise("load: " + msg +
                    ".\nError while trying to load behaviour '" + b +
                    "' in library '" + l + "'\n");
      }
    };
    //
    d.library = l;
    d.behaviour = b;
    d.function = fct;
    d.hypothesis = h;
    //
    if (lm.getMaterialKnowledgeType(l, b) != 1u) {
      raise("entry point '" + b + "' in library " + l + " is not a behaviour");
    }

    if (lm.getAPIVersion(l, b) != MGIS_BEHAVIOUR_API_VERSION) {
      std::string msg("unmatched API version\n");
      msg += "- the behaviour uses API version ";
      msg += std::to_string(lm.getAPIVersion(l, b)) + "\n";
      msg += "- mgis uses API version ";
      msg += std::to_string(MGIS_BEHAVIOUR_API_VERSION);
      raise(msg);
    }
    d.tfel_version = lm.getTFELVersion(l, b);
    d.unit_system = lm.getUnitSystem(l, b);
    d.source = lm.getSource(l, b);
    d.author = lm.getAuthor(l, b);
    d.date = lm.getDate(l, b);
    d.validator = lm.getValidator(l, b);
    d.build_id = lm.getBuildIdentifier(l, b);
    d.btype = [&l, &b, &lm, &raise] {
      /* - 0 : general behaviour
       * - 1 : strain based behaviour *
       * - 2 : standard finite strain behaviour *
       * - 3 : cohesive zone model */
      switch (lm.getBehaviourType(l, b)) {
        case 0:
          return BehaviourDescription::GENERALBEHAVIOUR;
        case 1:
          return BehaviourDescription::STANDARDSTRAINBASEDBEHAVIOUR;
        case 2:
          return BehaviourDescription::STANDARDFINITESTRAINBEHAVIOUR;
        case 3:
          return BehaviourDescription::COHESIVEZONEMODEL;
        default:
          break;
      }
      raise("unsupported behaviour type");
    }();
    // behaviour kinematic
    d.kinematic = [&l, &b, &lm, h, &raise] {
      /* - 0: undefined kinematic
       * - 1: standard small strain behaviour kinematic
       * - 2: cohesive zone model kinematic
       * - 3: standard finite strain kinematic (F-Cauchy)
       * - 4: ptest finite strain kinematic (eto-pk1)
       * - 5: Green-Lagrange strain
       * - 6: Miehe Apel Lambrecht logarithmic strain framework */
      switch (lm.getBehaviourKinematic(l, b)) {
        case 0:
          return BehaviourDescription::UNDEFINEDKINEMATIC;
        case 1:
          return BehaviourDescription::SMALLSTRAINKINEMATIC;
        case 2:
          return BehaviourDescription::COHESIVEZONEKINEMATIC;
        case 3:
          return BehaviourDescription::FINITESTRAINKINEMATIC_F_CAUCHY;
        case 4:
          if (((h != Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN) &&
               (h != Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS))) {
            raise(
                "invalid hypothesis for behaviour based on "
                "the eto-pk1 kinematic");
          }
          return BehaviourDescription::FINITESTRAINKINEMATIC_ETO_PK1;
      }
      raise("unsupported behaviour kinematic");
    }();
    // setting gradients and thermodynamic forces
    d.gradients = buildVariablesList(raise, lm.getGradientsNames(l, b, h),
                                     lm.getGradientsTypes(l, b, h));
    d.thermodynamic_forces =
        buildVariablesList(raise, lm.getThermodynamicForcesNames(l, b, h),
                           lm.getThermodynamicForcesTypes(l, b, h));
    raise_if(d.gradients.size() != d.thermodynamic_forces.size(),
             "the number of the gradients does not match "
             "the number of thermodynamic forces");
    switch (d.btype) {
      case BehaviourDescription::GENERALBEHAVIOUR:
        break;
      case BehaviourDescription::STANDARDSTRAINBASEDBEHAVIOUR:
        raise_if(d.kinematic != BehaviourDescription::SMALLSTRAINKINEMATIC,
                 "strain based behaviour must be associated with the "
                 "small strain kinematic hypothesis");
        checkGradientsAndThermodynamicForcesConsistency(
            raise, d.gradients, d.thermodynamic_forces,
            {"Strain", Variable::STENSOR, 1}, {"Stress", Variable::STENSOR, 1});
        break;
      case BehaviourDescription::COHESIVEZONEMODEL:
        if (d.kinematic != BehaviourDescription::COHESIVEZONEKINEMATIC) {
          raise("invalid kinematic assumption for cohesive zone model");
        }
        checkGradientsAndThermodynamicForcesConsistency(
            raise, d.gradients, d.thermodynamic_forces,
            {"OpeningDisplacement", Variable::VECTOR, 2},
            {"CohesiveForce", Variable::VECTOR, 2});
        break;
      case BehaviourDescription::STANDARDFINITESTRAINBEHAVIOUR:
        if (d.kinematic ==
            BehaviourDescription::FINITESTRAINKINEMATIC_F_CAUCHY) {
          checkGradientsAndThermodynamicForcesConsistency(
              raise, d.gradients, d.thermodynamic_forces,
              {"DeformationGradient", Variable::TENSOR, 3},
              {"Stress", Variable::STENSOR, 1});
        } else if (d.kinematic ==
                   BehaviourDescription::FINITESTRAINKINEMATIC_ETO_PK1) {
          if (((h != Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN) &&
               (h != Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS))) {
            raise(
                "invalid hypothesis for behaviour based on "
                "the eto-pk1 kinematic");
          }
          checkGradientsAndThermodynamicForcesConsistency(
              raise, d.gradients, d.thermodynamic_forces,
              {"Strain", Variable::STENSOR, 1},
              {"Stresss", Variable::STENSOR, 1});
        } else {
          raise(
              "invalid kinematic hypothesis for finite strain "
              "behaviour");
        }
        break;
      default:
        raise("unsupported behaviour type");
    }  // switch(d.btype)
    // behaviour symmetry
    d.symmetry = lm.getBehaviourSymmetry(l, b) == 0
                     ? BehaviourDescription::ISOTROPIC
                     : BehaviourDescription::ORTHOTROPIC;
    auto add_mp = [&d](const std::string &mp) {
      d.mps.push_back({mp, Variable::SCALAR});
    };
    if (lm.requiresStiffnessTensor(l, b, h)) {
      if (lm.getElasticStiffnessSymmetry(l, b) == 0) {
        add_mp("YoungModulus");
        add_mp("PoissonRatio");
      } else {
        if (d.symmetry != BehaviourDescription::ORTHOTROPIC) {
          raise(
              "load: the behaviour must be orthotropic "
              "for the elastic stiffness symmetry to be orthotropic");
        }
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
      if (d.symmetry == BehaviourDescription::ORTHOTROPIC) {
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
    d.isvs =
        buildVariablesList(raise, lm.getInternalStateVariablesNames(l, b, h),
                           lm.getInternalStateVariablesTypes(l, b, h));
    // external state variables
    if (lm.hasTemperatureBeenRemovedFromExternalStateVariables(l, b)) {
      d.esvs.push_back({"Temperature", Variable::SCALAR});
    }
    if (lm.hasExternalStateVariablesTypes(l, b, h)) {
      const auto esvs =
          buildVariablesList(raise, lm.getExternalStateVariablesNames(l, b, h),
                             lm.getExternalStateVariablesTypes(l, b, h));
      d.esvs.insert(d.esvs.end(), esvs.begin(), esvs.end());
    } else {
      // Prior to TFEL versions 3.4.4 and 4.1, external state variables were
      // only scalars
      for (const auto &esv : lm.getExternalStateVariablesNames(l, b, h)) {
        d.esvs.push_back({esv, Variable::SCALAR});
      }
    }
    // tangent operator blocks
    for (const auto &block : lm.getTangentOperatorBlocksNames(l, b, h)) {
      d.to_blocks.push_back(getJacobianBlockVariables(d, block));
    }
    d.computesStoredEnergy = lm.computesStoredEnergy(l, b, h);
    d.computesDissipatedEnergy = lm.computesDissipatedEnergy(l, b, h);
    //! parameters
    const auto pn = lm.getParametersNames(l, b, h);
    const auto pt = lm.getParametersTypes(l, b, h);
    raise_if(pn.size() != pt.size(),
             "inconsistent size between parameters' names and"
             "parameters' sizes");
    for (decltype(pn.size()) i = 0; i != pn.size(); ++i) {
      if (pt[i] == 0) {
        d.params.push_back(pn[i]);
      } else if (pt[i] == 1) {
        d.iparams.push_back(pn[i]);
      } else if (pt[i] == 2) {
        d.usparams.push_back(pn[i]);
      } else {
        raise("unsupported parameter type for parameter '" + pn[i] + "'");
      }
    }
  }  // end of loadBehaviourDescription

  bool isStandardFiniteStrainBehaviour(const std::string &l,
                                       const std::string &b) {
    auto &lm = mgis::LibrariesManager::get();
    return (lm.getBehaviourType(l, b) == 2) &&
           (lm.getBehaviourKinematic(l, b) == 3);
  }  // end of isStandardFiniteStrainBehaviour

  mgis::size_type getTangentOperatorArraySize(const BehaviourDescription &b) {
    auto s = mgis::size_type{};
    for (const auto &block : b.to_blocks) {
      s += getVariableSize(block.first, b.hypothesis) *
           getVariableSize(block.second, b.hypothesis);
    }
    return s;
  }  // end of getTangentOperatorArraySize

  template <>
  double getParameterDefaultValue<double>(const BehaviourDescription &b,
                                          const std::string &n) {
    auto &lm = mgis::LibrariesManager::get();
    return lm.getParameterDefaultValue(b.library, b.behaviour, b.hypothesis, n);
  }  // end of getParameterDefaultValue<double>

  template <>
  int getParameterDefaultValue<int>(const BehaviourDescription &b,
                                    const std::string &n) {
    auto &lm = mgis::LibrariesManager::get();
    return lm.getIntegerParameterDefaultValue(b.library, b.behaviour,
                                              b.hypothesis, n);
  }  // end of getParameterDefaultValue<int>

  template <>
  unsigned short getParameterDefaultValue<unsigned short>(
      const BehaviourDescription &b, const std::string &n) {
    auto &lm = mgis::LibrariesManager::get();
    return lm.getUnsignedShortParameterDefaultValue(b.library, b.behaviour,
                                                    b.hypothesis, n);
  }  // end of getParameterDefaultValue<unsigned short>

  bool hasBounds(const BehaviourDescription &b, const std::string &v) {
    auto &lm = mgis::LibrariesManager::get();
    return lm.hasBounds(b.library, b.behaviour, b.hypothesis, v);
  }  // end of hasBounds

  bool hasLowerBound(const BehaviourDescription &b, const std::string &v) {
    auto &lm = mgis::LibrariesManager::get();
    return lm.hasLowerBound(b.library, b.behaviour, b.hypothesis, v);
  }  // end of hasLowerBound

  bool hasUpperBound(const BehaviourDescription &b, const std::string &v) {
    auto &lm = mgis::LibrariesManager::get();
    return lm.hasUpperBound(b.library, b.behaviour, b.hypothesis, v);
  }  // end of hasUpperBound

  long double getLowerBound(const BehaviourDescription &b,
                            const std::string &v) {
    auto &lm = mgis::LibrariesManager::get();
    return lm.getLowerBound(b.library, b.behaviour, b.hypothesis, v);
  }  // end of getLowerBound

  long double getUpperBound(const BehaviourDescription &b,
                            const std::string &v) {
    auto &lm = mgis::LibrariesManager::get();
    return lm.getUpperBound(b.library, b.behaviour, b.hypothesis, v);
  }  // end of getUpperBound

  bool hasPhysicalBounds(const BehaviourDescription &b, const std::string &v) {
    auto &lm = mgis::LibrariesManager::get();
    return lm.hasPhysicalBounds(b.library, b.behaviour, b.hypothesis, v);
  }  // end of hasPhysicalBounds

  bool hasLowerPhysicalBound(const BehaviourDescription &b,
                             const std::string &v) {
    auto &lm = mgis::LibrariesManager::get();
    return lm.hasLowerPhysicalBound(b.library, b.behaviour, b.hypothesis, v);
  }  // end of hasLowerPhysicalBound

  bool hasUpperPhysicalBound(const BehaviourDescription &b,
                             const std::string &v) {
    auto &lm = mgis::LibrariesManager::get();
    return lm.hasUpperPhysicalBound(b.library, b.behaviour, b.hypothesis, v);
  }  // end of hasUpperPhysicalBound

  long double getLowerPhysicalBound(const BehaviourDescription &b,
                                    const std::string &v) {
    auto &lm = mgis::LibrariesManager::get();
    return lm.getLowerPhysicalBound(b.library, b.behaviour, b.hypothesis, v);
  }  // end of getLowerPhysicalBound

  long double getUpperPhysicalBound(const BehaviourDescription &b,
                                    const std::string &v) {
    auto &lm = mgis::LibrariesManager::get();
    return lm.getUpperPhysicalBound(b.library, b.behaviour, b.hypothesis, v);
  }  // end of getUpperPhysicalBound

  void print_markdown(std::ostream &,
                      const BehaviourDescription &,
                      const mgis::size_type) {}  // end of print_markdown

  std::vector<Variable> getBehaviourInitializeFunctionInputs(
      const std::string &l,
      const std::string &b,
      const std::string &i,
      const Hypothesis h) {
    auto &lm = mgis::LibrariesManager::get();
    auto raise = [&b, &l](const std::string &msg) {
      mgis::raise("load: " + msg + ".\nError while trying to load behaviour '" +
                  b + "' in library '" + l + "'\n");
    };
    return buildVariablesList(
        raise, lm.getBehaviourInitializeFunctionInputsNames(l, b, i, h),
        lm.getBehaviourInitializeFunctionInputsTypes(l, b, i, h));
  }  // end of getBehaviourInitializeFunctionInputs

  std::vector<Variable> getBehaviourPostProcessingOutputs(const std::string &l,
                                                          const std::string &b,
                                                          const std::string &i,
                                                          const Hypothesis h) {
    auto &lm = mgis::LibrariesManager::get();
    auto raise = [&b, &l](const std::string &msg) {
      mgis::raise("load: " + msg + ".\nError while trying to load behaviour '" +
                  b + "' in library '" + l + "'\n");
    };
    return buildVariablesList(
        raise, lm.getBehaviourPostProcessingOutputsNames(l, b, i, h),
        lm.getBehaviourPostProcessingOutputsTypes(l, b, i, h));
  }  // end of getBehaviourPostProcessingOutputs

}  // end of namespace mgis::behaviour
