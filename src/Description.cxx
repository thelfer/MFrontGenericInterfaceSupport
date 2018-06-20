/*!
 * \file   Description.cxx
 * \brief
 * \author HELFER Thomas 202608
 * \date   20 juin 2018
 */

#include "MFront/Raise.hxx"
#include "MFront/Behaviour/Description.hxx"

namespace mfront {

  namespace behaviour {

    Description::Description() = default;
    Description::Description(Description&&) = default;
    Description::Description(const Description&) = default;
    Description& Description::operator=(Description&&) = default;
    Description& Description::operator=(const Description&) = default;
    Description::~Description() = default;

    Description load(const std::string& l,
                     const std::string& b,
                     const Hypothesis h) {
      const auto fct = b + '_' + toString(h);
      auto d = Description{}; 

      d.library = l;
      d.behaviour = b;
      d.function = fct;
      d.hypothesis = h;

      //       d.tfel_version = elm.getTFELVersion(l, fct);
      //       d.source = elm.getSource(l, fct);
      //       d.btype = elm.getUMATBehaviourType(l, fct);
      //       d.kinematic = elm.getUMATBehaviourKinematic(l, fct);
      //       d.stype = elm.getUMATSymmetryType(l, fct);
      //       d.etype = elm.getUMATElasticSymmetryType(l, fct);
      //       d.isUPUIR =
      //           elm.isUMATBehaviourUsableInPurelyImplicitResolution(l, fct,
      //           h);
      //       d.mpnames = elm.getUMATMaterialPropertiesNames(l, fct, h);
      //       d.ivnames = elm.getUMATInternalStateVariablesNames(l, fct, h);
      //       d.ivtypes = elm.getUMATInternalStateVariablesTypes(l, fct, h);
      //       d.evnames = elm.getUMATExternalStateVariablesNames(l, fct, h);
      //       //! parameters
      //       const auto pn = elm.getUMATParametersNames(l, fct, h);
      //       const auto pt = elm.getUMATParametersTypes(l, fct, h);
      //       throw_if(
      //           pn.size() != pt.size(),
      //           "inconsistent size between parameters' names and parameters'
      //           sizes");
      //       for (decltype(pn.size()) i = 0; i != pn.size(); ++i) {
      //         if (pt[i] == 0) {
      //           d.pnames.push_back(pn[i]);
      //         } else if (pt[i] == 1) {
      //           d.ipnames.push_back(pn[i]);
      //         } else if (pt[i] == 2) {
      //           d.upnames.push_back(pn[i]);
      //         } else {
      //           throw_if(true,
      //                    "unsupported parameter type for parameter '" + pn[i]
      //                    + "'");
      //         }
      //       }
      //       //! additional parameters
      //       d.requiresStiffnessTensor =
      //           elm.getUMATRequiresStiffnessTensor(l, fct, h);
      //       d.requiresThermalExpansionCoefficientTensor =
      //           elm.getUMATRequiresThermalExpansionCoefficientTensor(l, fct,
      //           h);

      return d;
    }  // end of load

  }  // end of namespace behaviour

}  // end of namespace mfront
