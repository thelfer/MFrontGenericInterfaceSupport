/*!
 * \file   Behaviour.cxx
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

#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/LibrariesManager.hxx"
#include "MGIS/Raise.hxx"

namespace mgis {

  namespace behaviour {

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
        switch (types[i]) {
          case 0:
            vars.push_back({names[i], Variable::SCALAR});
            break;
          case 1:
            vars.push_back({names[i], Variable::STENSOR});
            break;
          case 2:
            vars.push_back({names[i], Variable::VECTOR});
            break;
          case 3:
            vars.push_back({names[i], Variable::TENSOR});
            break;
          default:
            raise("unsupported internal state variable type");
        }
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

    Behaviour::Behaviour() = default;
    Behaviour::Behaviour(Behaviour &&) = default;
    Behaviour::Behaviour(const Behaviour &) = default;
    Behaviour &Behaviour::operator=(Behaviour &&) = default;
    Behaviour &Behaviour::operator=(const Behaviour &) = default;
    Behaviour::~Behaviour() = default;

    Behaviour load(const std::string &l,
                   const std::string &b,
                   const Hypothesis h) {
      auto &lm = mgis::LibrariesManager::get();
      const auto fct = b + '_' + toString(h);
      auto raise = [&b, &l](const std::string &msg) {
        mgis::raise("load: " + msg +
                    ".\nError while trying to load behaviour '" + b +
                    "' in library '" + l + "'\n");
      };
      auto raise_if = [&b, &l](const bool c, const std::string &msg) {
        if (c) {
          mgis::raise("load: " + msg +
                      ".\nError while trying to load behaviour '" + b +
                      "' in library '" + l + "'\n");
        }
      };
      auto d = Behaviour{};

      d.library = l;
      d.behaviour = b;
      d.function = fct;
      d.hypothesis = h;
      d.b = lm.getBehaviour(l, b, h);

      if (lm.getMaterialKnowledgeType(l, b) != 1u) {
        raise("entry point '" + b + "' in library " + l +
              " is not a behaviour");
      }

      d.tfel_version = lm.getTFELVersion(l, b);
      d.source = lm.getSource(l, b);
      d.btype = [&l, &b, &lm, &raise] {
        /* - 0 : general behaviour
         * - 1 : strain based behaviour *
         * - 2 : standard finite strain behaviour *
         * - 3 : cohesive zone model */
        switch (lm.getBehaviourType(l, b)) {
          case 0:
            return Behaviour::GENERALBEHAVIOUR;
            break;
          case 1:
            return Behaviour::STANDARDSTRAINBASEDBEHAVIOUR;
            break;
          case 2:
            return Behaviour::STANDARDFINITESTRAINBEHAVIOUR;
            break;
          case 3:
            return Behaviour::COHESIVEZONEMODEL;
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
            return Behaviour::UNDEFINEDKINEMATIC;
            break;
          case 1:
            return Behaviour::SMALLSTRAINKINEMATIC;
            break;
          case 2:
            return Behaviour::COHESIVEZONEKINEMATIC;
            break;
          case 3:
            return Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY;
            break;
          case 4:
            if (((h != Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN) &&
                 (h != Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS))) {
              raise(
                  "invalid hypothesis for behaviour based on "
                  "the eto-pk1 kinematic");
            }
            return Behaviour::FINITESTRAINKINEMATIC_ETO_PK1;
            break;
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
        case Behaviour::GENERALBEHAVIOUR:
          break;
        case Behaviour::STANDARDSTRAINBASEDBEHAVIOUR:
          raise_if(d.kinematic != Behaviour::SMALLSTRAINKINEMATIC,
                   "strain based behaviour must be associated with the "
                   "small strain kinematic hypothesis");
          checkGradientsAndThermodynamicForcesConsistency(
              raise, d.gradients, d.thermodynamic_forces,
              {"Strain", Variable::STENSOR}, {"Stress", Variable::STENSOR});
          break;
        case Behaviour::COHESIVEZONEMODEL:
          if (d.kinematic != Behaviour::COHESIVEZONEKINEMATIC) {
            raise("invalid kinematic assumption for cohesive zone model");
          }
          checkGradientsAndThermodynamicForcesConsistency(
              raise, d.gradients, d.thermodynamic_forces,
              {"OpeningDisplacement", Variable::VECTOR},
              {"CohesiveForce", Variable::VECTOR});
          break;
        case Behaviour::STANDARDFINITESTRAINBEHAVIOUR:
          if (d.kinematic == Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY) {
            checkGradientsAndThermodynamicForcesConsistency(
                raise, d.gradients, d.thermodynamic_forces,
                {"DeformationGradient", Variable::TENSOR},
                {"Stress", Variable::STENSOR});
          } else if (d.kinematic == Behaviour::FINITESTRAINKINEMATIC_ETO_PK1) {
            if (((h != Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN) &&
                 (h != Hypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS))) {
              raise(
                  "invalid hypothesis for behaviour based on "
                  "the eto-pk1 kinematic");
            }
            checkGradientsAndThermodynamicForcesConsistency(
                raise, d.gradients, d.thermodynamic_forces,
                {"Strain", Variable::STENSOR}, {"Stresss", Variable::STENSOR});
          } else {
            raise(
                "invalid kinematic hypothesis for finite strain "
                "behaviour");
          }
          break;
        default:
          raise("unsupported behaviour type");
      };
      // behaviour symmetry
      d.symmetry = lm.getBehaviourSymmetry(l, b) == 0 ? Behaviour::ISOTROPIC
                                                      : Behaviour::ORTHOTROPIC;
      auto add_mp = [&d](const std::string &mp) {
        d.mps.push_back({mp, Variable::SCALAR});
      };
      if (lm.requiresStiffnessTensor(l, b, h)) {
        if (lm.getElasticStiffnessSymmetry(l, b) == 0) {
          add_mp("YoungModulus");
          add_mp("PoissonRatio");
        } else {
          if (d.symmetry != Behaviour::ORTHOTROPIC) {
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
        if (d.symmetry == Behaviour::ORTHOTROPIC) {
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
      d.esvs.push_back({"Temperature", Variable::SCALAR});
      for (const auto &esv : lm.getExternalStateVariablesNames(l, b, h)) {
        d.esvs.push_back({esv, Variable::SCALAR});
      }
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
      return d;
    }  // end of load

    void setParameter(const Behaviour &b,
                      const std::string &n,
                      const double v) {
      auto &lm = mgis::LibrariesManager::get();
      lm.setParameter(b.library, b.function, b.hypothesis, n, v);
    } // end of setParameter

    void setParameter(const Behaviour &b,
                      const std::string &n,
                      const int v) {
      auto &lm = mgis::LibrariesManager::get();
      lm.setParameter(b.library, b.function, b.hypothesis, n, v);
    } // end of setParameter

    void setParameter(const Behaviour &b,
                      const std::string &n,
                      const unsigned short v) {
      auto &lm = mgis::LibrariesManager::get();
      lm.setParameter(b.library, b.function, b.hypothesis, n, v);
    } // end of setParameter

    template <>
    double getParameterDefaultValue<double>(const Behaviour &b,
                                            const std::string &n) {
      auto &lm = mgis::LibrariesManager::get();
      return lm.getParameterDefaultValue(b.library, b.function, b.hypothesis,
                                         n);
    }  // end of getParameterDefaultValue<double>

    template <>
    int getParameterDefaultValue<int>(const Behaviour &b,
                                      const std::string &n) {
      auto &lm = mgis::LibrariesManager::get();
      return lm.getIntegerParameterDefaultValue(b.library, b.function,
                                                b.hypothesis, n);
    }  // end of getParameterDefaultValue<int>

    template <>
    unsigned short getParameterDefaultValue<unsigned short>(
        const Behaviour &b, const std::string &n) {
      auto &lm = mgis::LibrariesManager::get();
      return lm.getUnsignedShortParameterDefaultValue(b.library, b.function,
                                                      b.hypothesis, n);
    }  // end of getParameterDefaultValue<unsigned short>

    bool hasBounds(const Behaviour &b, const std::string &v) {
      auto &lm = mgis::LibrariesManager::get();
      return lm.hasBounds(b.library, b.function, b.hypothesis, v);
    }  // end of hasBounds

    bool hasLowerBound(const Behaviour &b, const std::string &v) {
      auto &lm = mgis::LibrariesManager::get();
      return lm.hasLowerBound(b.library, b.function, b.hypothesis, v);
    }  // end of hasLowerBound

    bool hasUpperBound(const Behaviour &b, const std::string &v) {
      auto &lm = mgis::LibrariesManager::get();
      return lm.hasUpperBound(b.library, b.function, b.hypothesis, v);
    }  // end of hasUpperBound

    long double getLowerBound(const Behaviour &b, const std::string &v) {
      auto &lm = mgis::LibrariesManager::get();
      return lm.getLowerBound(b.library, b.function, b.hypothesis, v);
    }  // end of getLowerBound

    long double getUpperBound(const Behaviour &b, const std::string &v) {
      auto &lm = mgis::LibrariesManager::get();
      return lm.getUpperBound(b.library, b.function, b.hypothesis, v);
    }  // end of getUpperBound

    bool hasPhysicalBounds(const Behaviour &b, const std::string &v) {
      auto &lm = mgis::LibrariesManager::get();
      return lm.hasPhysicalBounds(b.library, b.function, b.hypothesis, v);
    }  // end of hasPhysicalBounds

    bool hasLowerPhysicalBound(const Behaviour &b, const std::string &v) {
      auto &lm = mgis::LibrariesManager::get();
      return lm.hasLowerPhysicalBound(b.library, b.function, b.hypothesis, v);
    }  // end of hasLowerPhysicalBound

    bool hasUpperPhysicalBound(const Behaviour &b, const std::string &v) {
      auto &lm = mgis::LibrariesManager::get();
      return lm.hasUpperPhysicalBound(b.library, b.function, b.hypothesis, v);
    }  // end of hasUpperPhysicalBound

    long double getLowerPhysicalBound(const Behaviour &b,
                                      const std::string &v) {
      auto &lm = mgis::LibrariesManager::get();
      return lm.getLowerPhysicalBound(b.library, b.function, b.hypothesis, v);
    }  // end of getLowerPhysicalBound

    long double getUpperPhysicalBound(const Behaviour &b,
                                      const std::string &v) {
      auto &lm = mgis::LibrariesManager::get();
      return lm.getUpperPhysicalBound(b.library, b.function, b.hypothesis, v);
    }  // end of getUpperPhysicalBound

  }  // end of namespace behaviour

}  // end of namespace mgis
