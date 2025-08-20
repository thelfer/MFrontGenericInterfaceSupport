/*!
 * \file   LibrariesManager.cxx
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

#include <cctype>
#include <cstring>
#include <iterator>
#include <algorithm>
#if !((defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__))
#include <dlfcn.h>
#endif /* !((defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)) */

#include "MGIS/LibrariesManager.hxx"
#include "MGIS/Raise.hxx"

namespace mgis {

#if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
  // code retrieved from
  // http://www.codeproject.com/Tips/479880/GetLastError-as-std-string
  static std::string getLastWin32Error() {
    const DWORD error = GetLastError();
    if (error) {
      LPVOID lpMsgBuf;
      DWORD bufLen = FormatMessage(
          FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
              FORMAT_MESSAGE_IGNORE_INSERTS,
          nullptr, error, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
          (LPTSTR)&lpMsgBuf, 0, nullptr);
      if (bufLen) {
        LPCSTR lpMsgStr = (LPTSTR)lpMsgBuf;
        std::string result(lpMsgStr, lpMsgStr + bufLen);
        LocalFree(lpMsgBuf);
        return result;
      }
    }
    return std::string();
  }
#endif /*  (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */

  static std::string getErrorMessage() {
#if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
    return getLastWin32Error();
#else
    const auto e = ::dlerror();
    if (e != nullptr) {
      return std::string(e);
    }
    return "";
#endif /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */
  }    // end of getErrorMessage

  static LibrariesManager::libhandler load_library(const std::string &l) {
#if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
    return ::LoadLibrary(TEXT(l.c_str()));
#else
    return ::dlopen(l.c_str(), RTLD_NOW);
#endif /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */
  }    // end of load_library

  static std::pair<LibrariesManager::libhandler, std::string> try_open(
      const std::string &l) {
    auto starts_with = [](const std::string &s1, const char *const s2) {
      const auto ls2 = std::strlen(s2);
      return ((s1.size() >= ls2) && (std::equal(s2, s2 + ls2, s1.begin())));
    };  // end of starts_with
    auto ends_with = [](const std::string &s1, const char *const s2) {
      const auto ls2 = std::strlen(s2);
      if (!(s1.size() >= ls2)) {
        return false;
      }
      return std::equal(s2, s2 + ls2, s1.begin() + (s1.size() - ls2));
    };  // end of ends_with
#if (defined(macintosh) || defined(Macintosh) || \
     (defined(__APPLE__) && defined(__MACH__)))
    const char *const ext = ".dylib";
#elif (defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__))
    const char *const ext = ".dll";
#else
    const char *const ext = ".so";
#endif
    auto ln = l;
    auto lib = load_library(l);
#if defined(__CYGWIN__)
    if ((lib == nullptr) && (!starts_with(l, "cyg"))) {
      ln = "cyg" + l;
      lib = load_library(ln);
      if (lib == nullptr) {
        if (!ends_with(l, ext)) {
          ln = "cyg" + l + ext;
          lib = load_library(ln);
        }
      }
    }
#endif
#if !(defined(_WIN32) || defined(_WIN64))
    if ((lib == nullptr) && (!starts_with(l, "lib"))) {
      ln = "lib" + l;
      lib = load_library(ln);
      if (lib == nullptr) {
        if (!ends_with(l, ext)) {
          ln = "lib" + l + ext;
          lib = load_library(ln);
        }
      }
    }
#endif
    if ((lib == nullptr) && (!ends_with(l, ext))) {
      ln = l + ext;
      lib = load_library(ln);
    }
    // retrieving the initial error message
    if (lib == nullptr) {
      ln = l;
      lib = load_library(ln);
    }
    return {lib, ln};
  }  // end of try_open

  static std::string decomposeVariableName1(const std::string &n) {
    auto throw_if = [](const bool c, const std::string &m) {
      mgis::raise_if(c, "mgis::decomposeVariableName: " + m);
    };
    auto p = n.cbegin();
    auto pe = n.cend();
    while ((p != pe) && (*p != '[')) {
      ++p;
    }
    if (p == pe) {
      return n;
    }
    auto r = std::string{n.cbegin(), p};
    ++p;
    throw_if(p == pe, "unexpected end of string 'n'");
    throw_if(!std::isdigit(*p), "unexpected a digit 'n'");
    r += "__";
    while ((p != pe) && (std::isdigit(*p))) {
      r.push_back(*p);
      ++p;
    }
    throw_if(p == pe, "unexpected end of string '" + n + "'");
    throw_if(*p != ']', "invalid variable name '" + n + "'");
    ++p;
    throw_if(p != pe, "invalid variable name '" + n + "'");
    r += "__";
    return r;
  }  // end of decomposeVariableName1

  static std::string decomposeVariableName2(const std::string &n) {
    auto throw_if = [](const bool c, const std::string &m) {
      mgis::raise_if(c, "mgis::decomposeVariableName: " + m);
    };
    auto p = n.cbegin();
    auto pe = n.cend();
    while ((p != pe) && (*p != '[')) {
      ++p;
    }
    if (p == pe) {
      return n;
    }
    auto r = std::string{n.cbegin(), p};
    ++p;
    throw_if(p == pe, "unexpected end of string 'n'");
    throw_if(!std::isdigit(*p), "unexpected a digit 'n'");
    r += "_mfront_index_";
    while ((p != pe) && (std::isdigit(*p))) {
      r.push_back(*p);
      ++p;
    }
    throw_if(p == pe, "unexpected end of string '" + n + "'");
    throw_if(*p != ']', "invalid variable name '" + n + "'");
    ++p;
    throw_if(p != pe, "invalid variable name '" + n + "'");
    return r;
  }  // end of decomposeVariableName2

  static std::pair<std::string, std::string> decomposeVariableName(
      const std::string &n) {
    return {decomposeVariableName1(n), decomposeVariableName2(n)};
  }  // end of decomposeVariableName

  LibrariesManager &LibrariesManager::get() {
    static LibrariesManager lm;
    return lm;
  }  // end of get

  LibrariesManager::LibrariesManager() = default;

  mgis::material_property::MaterialPropertyFctPtr
  LibrariesManager::getMaterialProperty(const std::string &l,
                                        const std::string &mp) {
    const auto p = this->getSymbolAddress(l, mp);
    if (p == nullptr) {
      mgis::raise(
          "LibrariesManager::getMaterialProperty: can't load material property "
          "'" +
          mp + "' in library '" + l + "'");
    }
    return reinterpret_cast<mgis::material_property::MaterialPropertyFctPtr>(p);
  }  // end of getMaterialProperty

  std::vector<std::string> LibrariesManager::getMaterialPropertyInputsNames(
      const std::string &l, const std::string &mp) {
    return this->getNames(l, mp, "args");
  }  // end of getMaterialPropertyInputsNames

  std::string LibrariesManager::getMaterialPropertyOutputName(
      const std::string &l, const std::string &mp) {
    return *(this->extract<const char *const>(l, mp + "_output"));
  }  // end of getMaterialPropertyOutputName

  std::vector<std::string> LibrariesManager::getBehaviourInitializeFunctions(
      const std::string &l, const std::string &b, const Hypothesis h) {
    return this->getNames(l, b, h, "InitializeFunctions");
  }  // end of getBehaviourInitializeFunctions

  mgis::behaviour::BehaviourInitializeFctPtr
  LibrariesManager::getBehaviourInitializeFunction(const std::string &l,
                                                   const std::string &b,
                                                   const std::string &f,
                                                   const Hypothesis h) {
    const auto hn = toString(h);
    const auto p =
        this->getSymbolAddress(l, b + "_" + hn + "_InitializeFunction_" + f);
    if (p == nullptr) {
      mgis::raise(
          "LibrariesManager::getBehaviourInitializeFunction: "
          "can't load initialize function '" +
          f + "' for behaviour '" + b + "' in library '" + l +
          "' for hypothesis '" + hn + "'");
    }
    return reinterpret_cast<mgis::behaviour::BehaviourInitializeFctPtr>(p);
  }  // end of getBehaviour

  std::vector<std::string>
  LibrariesManager::getBehaviourInitializeFunctionInputsNames(
      const std::string &l,
      const std::string &b,
      const std::string &p,
      const Hypothesis h) {
    return this->getNames(l, b, h, "InitializeFunction_" + p + "_Inputs");
  }  // end of getBehaviourInitializeFunctionInputsNames

  std::vector<int> LibrariesManager::getBehaviourInitializeFunctionInputsTypes(
      const std::string &l,
      const std::string &b,
      const std::string &p,
      const Hypothesis h) {
    auto types = std::vector<int>{};
    const auto hn = toString(h);
    const auto outputs = "InitializeFunction_" + p + "_Inputs";
    const auto outputs_types = "InitializeFunction_" + p + "_InputsTypes";
    const auto nb = *(this->extract<unsigned short>(
        l, b + "_" + hn + "_n" + outputs, b + "_n" + outputs));
    const auto res = this->extract<const int>(
        l, b + "_" + hn + '_' + outputs_types, b + '_' + outputs_types);
    std::copy(res, res + nb, std::back_inserter(types));
    return types;
  }  // end of getBehaviourInitializeFunctionInputsTypes

  mgis::behaviour::BehaviourFctPtr LibrariesManager::getBehaviour(
      const std::string &l, const std::string &b, const Hypothesis h) {
    const auto hn = toString(h);
    const auto p = this->getSymbolAddress(l, b + "_" + hn);
    if (p == nullptr) {
      mgis::raise(
          "LibrariesManager::getBehaviour: "
          "can't load behaviour '" +
          b + "' in library '" + l + "' for hypothesis '" + hn + "'");
    }
    return reinterpret_cast<mgis::behaviour::BehaviourFctPtr>(p);
  }  // end of getBehaviour

  std::vector<std::string> LibrariesManager::getBehaviourPostProcessings(
      const std::string &l, const std::string &b, const Hypothesis h) {
    return this->getNames(l, b, h, "PostProcessings");
  }  // end of getBehaviourPostProcessings

  mgis::behaviour::BehaviourPostProcessingFctPtr
  LibrariesManager::getBehaviourPostProcessing(const std::string &l,
                                               const std::string &b,
                                               const std::string &f,
                                               const Hypothesis h) {
    const auto hn = toString(h);
    const auto p =
        this->getSymbolAddress(l, b + "_" + hn + "_PostProcessing_" + f);
    if (p == nullptr) {
      mgis::raise(
          "LibrariesManager::getBehaviourPostProcessing: "
          "can't load post-processing '" +
          f + "' for behaviour '" + b + "' in library '" + l +
          "' for hypothesis '" + hn + "'");
    }
    return reinterpret_cast<mgis::behaviour::BehaviourPostProcessingFctPtr>(p);
  }  // end of getBehaviourPostProcessing

  std::vector<std::string>
  LibrariesManager::getBehaviourPostProcessingOutputsNames(const std::string &l,
                                                           const std::string &b,
                                                           const std::string &p,
                                                           const Hypothesis h) {
    return this->getNames(l, b, h, "PostProcessing_" + p + "_Outputs");
  }  // end of getBehaviourPostProcessingOutputsNames

  std::vector<int> LibrariesManager::getBehaviourPostProcessingOutputsTypes(
      const std::string &l,
      const std::string &b,
      const std::string &p,
      const Hypothesis h) {
    auto types = std::vector<int>{};
    const auto hn = toString(h);
    const auto outputs = "PostProcessing_" + p + "_Outputs";
    const auto outputs_types = "PostProcessing_" + p + "_OutputsTypes";
    const auto nb = *(this->extract<unsigned short>(
        l, b + "_" + hn + "_n" + outputs, b + "_n" + outputs));
    const auto res = this->extract<const int>(
        l, b + "_" + hn + '_' + outputs_types, b + '_' + outputs_types);
    std::copy(res, res + nb, std::back_inserter(types));
    return types;
  }  // end of getBehaviourPostProcessingOutputsTypes

  mgis::behaviour::RotateBehaviourGradientsFctPtr
  LibrariesManager::getRotateBehaviourGradientsFunction(const std::string &l,
                                                        const std::string &b,
                                                        const Hypothesis h) {
    const auto hn = toString(h);
    const auto f = b + "_" + hn + "_rotateGradients";
    const auto p = this->getSymbolAddress(l, f);
    if (p == nullptr) {
      mgis::raise(
          "LibrariesManager::getRotateBehaviourGradientsFunction: "
          "can't load gradients' rotation function '" +
          f +
          "' "
          "for behaviour '" +
          b + "' in library '" + l + "' for hypothesis '" + hn + "'");
    }
    return reinterpret_cast<mgis::behaviour::RotateBehaviourGradientsFctPtr>(p);
  }  // end of getRotateBehaviourGradientsFunction

  mgis::behaviour::RotateArrayOfBehaviourGradientsFctPtr
  LibrariesManager::getRotateArrayOfBehaviourGradientsFunction(
      const std::string &l, const std::string &b, const Hypothesis h) {
    const auto hn = toString(h);

    const auto f = b + "_" + hn + "_rotateArrayOfGradients";
    const auto p = this->getSymbolAddress(l, f);
    if (p == nullptr) {
      mgis::raise(
          "LibrariesManager::getRotateArrayOfBehaviourGradientsFunction: "
          "can't load gradients' rotation function '" +
          f +
          "' "
          "for behaviour '" +
          b + "' in library '" + l + "' for hypothesis '" + hn + "'");
    }
    return reinterpret_cast<
        mgis::behaviour::RotateArrayOfBehaviourGradientsFctPtr>(p);
  }  // end of getRotateArrayOfBehaviourGradientsFunction

  mgis::behaviour::RotateBehaviourThermodynamicForcesFctPtr
  LibrariesManager::getRotateBehaviourThermodynamicForcesFunction(
      const std::string &l, const std::string &b, const Hypothesis h) {
    const auto hn = toString(h);
    const auto f = b + "_" + hn + "_rotateThermodynamicForces";
    const auto p = this->getSymbolAddress(l, f);
    if (p == nullptr) {
      mgis::raise(
          "LibrariesManager::getRotateBehaviourThermodynamicForcesFunction: "
          "can't load thermodynamic forces' rotation function '" +
          f + "' for behaviour '" + b + "' in library '" + l +
          "' for hypothesis '" + hn + "'");
    }
    return reinterpret_cast<
        mgis::behaviour::RotateBehaviourThermodynamicForcesFctPtr>(p);
  }  // end of getRotateBehaviourThermodynamicForcesFunction

  mgis::behaviour::RotateArrayOfBehaviourThermodynamicForcesFctPtr
  LibrariesManager::getRotateArrayOfBehaviourThermodynamicForcesFunction(
      const std::string &l, const std::string &b, const Hypothesis h) {
    const auto hn = toString(h);
    const auto f = b + "_" + hn + "_rotateArrayOfThermodynamicForces";
    const auto p = this->getSymbolAddress(l, f);
    if (p == nullptr) {
      mgis::raise(
          "LibrariesManager::"
          "getRotateArrayOfBehaviourThermodynamicForcesFunction: "
          "can't load array of thermodynamic forces' rotation function '" +
          f + "' for behaviour '" + b + "' in library '" + l +
          "' for hypothesis '" + hn + "'");
    }
    return reinterpret_cast<
        mgis::behaviour::RotateArrayOfBehaviourThermodynamicForcesFctPtr>(p);
  }  // end of getRotateArrayOfBehaviourThermodynamicForcesFunction

  mgis::behaviour::RotateBehaviourThermodynamicForcesFctPtr
  LibrariesManager::getRotateBehaviourThermodynamicForcesFunction(
      const std::string &l,
      const std::string &b,
      const Hypothesis h,
      const mgis::behaviour::FiniteStrainBehaviourOptions::StressMeasure s) {
    const auto hn = toString(h);
    const auto suffix = [&s]() -> std::string {
      if (s == mgis::behaviour::FiniteStrainBehaviourOptions::CAUCHY) {
        return "CauchyStress";
      } else if (s == mgis::behaviour::FiniteStrainBehaviourOptions::PK2) {
        return "PK2Stress";
      } else if (s != mgis::behaviour::FiniteStrainBehaviourOptions::PK1) {
        mgis::raise(
            "LibrariesManager::getRotateBehaviourTangentOperatorBlocksFunction:"
            " unsupported stress measure");
      }
      return "PK1Stress";
    }();
    const auto f = b + "_" + hn + "_rotateThermodynamicForces_" + suffix;
    const auto p = this->getSymbolAddress(l, f);
    if (p == nullptr) {
      mgis::raise(
          "LibrariesManager::getRotateBehaviourThermodynamicForcesFunction: "
          "can't load gradients' rotation function '" +
          f + "' for behaviour '" + b + "' in library '" + l +
          "' for hypothesis '" + hn + "'");
    }
    return reinterpret_cast<
        mgis::behaviour::RotateBehaviourThermodynamicForcesFctPtr>(p);
  }  // end of getRotateBehaviourThermodynamicForcesFunction

  mgis::behaviour::RotateArrayOfBehaviourThermodynamicForcesFctPtr
  LibrariesManager::getRotateArrayOfBehaviourThermodynamicForcesFunction(
      const std::string &l,
      const std::string &b,
      const Hypothesis h,
      const mgis::behaviour::FiniteStrainBehaviourOptions::StressMeasure s) {
    const auto hn = toString(h);
    const auto suffix = [&s]() -> std::string {
      if (s == mgis::behaviour::FiniteStrainBehaviourOptions::CAUCHY) {
        return "CauchyStress";
      } else if (s == mgis::behaviour::FiniteStrainBehaviourOptions::PK2) {
        return "PK2Stress";
      }
      mgis::raise_if(s != mgis::behaviour::FiniteStrainBehaviourOptions::PK1,
                     "LibrariesManager::"
                     "getRotateArrayOfBehaviourTangentOperatorBlocksFunction:"
                     " unsupported stress measure");
      return "PK1Stress";
    }();
    const auto f = b + "_" + hn + "_rotateArrayOfThermodynamicForces_" + suffix;
    const auto p = this->getSymbolAddress(l, f);
    if (p == nullptr) {
      mgis::raise(
          "LibrariesManager::"
          "getRotateArrayOfBehaviourThermodynamicForcesFunction: "
          "can't load thermodynamic forces' rotation function '" +
          f + "' for behaviour '" + b + "' in library '" + l +
          "' for hypothesis '" + hn + "'");
    }
    return reinterpret_cast<
        mgis::behaviour::RotateArrayOfBehaviourThermodynamicForcesFctPtr>(p);
  }  // end of getRotateArrayOfBehaviourThermodynamicForcesFunction

  mgis::behaviour::RotateBehaviourTangentOperatorBlocksFctPtr
  LibrariesManager::getRotateBehaviourTangentOperatorBlocksFunction(
      const std::string &l, const std::string &b, const Hypothesis h) {
    const auto hn = toString(h);
    const auto f = b + "_" + hn + "_rotateTangentOperatorBlocks";
    const auto p = this->getSymbolAddress(l, f);
    if (p == nullptr) {
      mgis::raise(
          "LibrariesManager::getRotateBehaviourTangentOperatorBlocksFunction: "
          "can't load tangent operator blocks' rotation function '" +
          f + "' for behaviour '" + b + "' in library '" + l +
          "' for hypothesis '" + hn + "'");
    }
    return reinterpret_cast<
        mgis::behaviour::RotateBehaviourTangentOperatorBlocksFctPtr>(p);
  }  // end of getRotateBehaviourTangentOperatorBlocksFunction

  mgis::behaviour::RotateArrayOfBehaviourTangentOperatorBlocksFctPtr
  LibrariesManager::getRotateArrayOfBehaviourTangentOperatorBlocksFunction(
      const std::string &l, const std::string &b, const Hypothesis h) {
    const auto hn = toString(h);
    const auto f = b + "_" + hn + "_rotateArrayOfTangentOperatorBlocks";
    const auto p = this->getSymbolAddress(l, f);
    if (p == nullptr) {
      mgis::raise(
          "LibrariesManager::"
          "getRotateArrayOfBehaviourTangentOperatorBlocksFunction: "
          "can't load tangent operator blocks' rotation function '" +
          f + "' for behaviour '" + b + "' in library '" + l +
          "' for hypothesis '" + hn + "'");
    }
    return reinterpret_cast<
        mgis::behaviour::RotateArrayOfBehaviourTangentOperatorBlocksFctPtr>(p);
  }  // end of getRotateArrayOfBehaviourTangentOperatorBlocksFunction

  mgis::behaviour::RotateBehaviourTangentOperatorBlocksFctPtr
  LibrariesManager::getRotateBehaviourTangentOperatorBlocksFunction(
      const std::string &l,
      const std::string &b,
      const Hypothesis h,
      const mgis::behaviour::FiniteStrainBehaviourOptions::TangentOperator t) {
    const auto hn = toString(h);
    const auto suffix = [&t]() -> std::string {
      if (t == mgis::behaviour::FiniteStrainBehaviourOptions::DSIG_DF) {
        return "dsig_dF";
      } else if (t == mgis::behaviour::FiniteStrainBehaviourOptions::DS_DEGL) {
        return "dPK2_dEGL";
      } else if (t == mgis::behaviour::FiniteStrainBehaviourOptions::DTAU_DDF) {
        return "dtau_ddF";
      } else if (t != mgis::behaviour::FiniteStrainBehaviourOptions::DPK1_DF) {
        mgis::raise(
            "LibrariesManager::getRotateBehaviourTangentOperatorBlocksFunction:"
            " unsupported tangent operator type");
      }
      return "dPK1_dF";
    }();
    const auto f = b + "_" + hn + "_rotateTangentOperatorBlocks_" + suffix;
    const auto p = this->getSymbolAddress(l, f);
    if (p == nullptr) {
      mgis::raise(
          "LibrariesManager::getRotateBehaviourTangentOperatorBlocksFunction: "
          "can't load tangent operator blocks' rotation function '" +
          f + "' for behaviour '" + b + "' in library '" + l +
          "' for hypothesis '" + hn + "'");
    }
    return reinterpret_cast<
        mgis::behaviour::RotateBehaviourTangentOperatorBlocksFctPtr>(p);
  }  // end of getRotateBehaviourTangentOperatorBlocksFunction

  mgis::behaviour::RotateArrayOfBehaviourTangentOperatorBlocksFctPtr
  LibrariesManager::getRotateArrayOfBehaviourTangentOperatorBlocksFunction(
      const std::string &l,
      const std::string &b,
      const Hypothesis h,
      const mgis::behaviour::FiniteStrainBehaviourOptions::TangentOperator t) {
    const auto hn = toString(h);
    const auto suffix = [&t]() -> std::string {
      if (t == mgis::behaviour::FiniteStrainBehaviourOptions::DSIG_DF) {
        return "dsig_dF";
      } else if (t == mgis::behaviour::FiniteStrainBehaviourOptions::DS_DEGL) {
        return "dPK2_dEGL";
      } else if (t == mgis::behaviour::FiniteStrainBehaviourOptions::DTAU_DDF) {
        return "dtau_ddF";
      } else if (t != mgis::behaviour::FiniteStrainBehaviourOptions::DPK1_DF) {
        mgis::raise(
            "LibrariesManager::"
            "getRotateArrayOfBehaviourTangentOperatorBlocksFunction:"
            " unsupported tangent operator type");
      }
      return "dPK1_dF";
    }();
    const auto f =
        b + "_" + hn + "_rotateArrayOfTangentOperatorBlocks_" + suffix;
    const auto p = this->getSymbolAddress(l, f);
    if (p == nullptr) {
      mgis::raise(
          "LibrariesManager::"
          "getRotateArrayOfBehaviourTangentOperatorBlocksFunction: "
          "can't load tangent operator blocks' rotation function '" +
          f + "' for behaviour '" + b + "' in library '" + l +
          "' for hypothesis '" + hn + "'");
    }
    return reinterpret_cast<
        mgis::behaviour::RotateArrayOfBehaviourTangentOperatorBlocksFctPtr>(p);
  }  // end of getRotateArrayOfBehaviourTangentOperatorBlocksFunction

  std::string LibrariesManager::getTFELVersion(const std::string &l,
                                               const std::string &n) {
    const auto p = this->getSymbolAddress(l, n + "_tfel_version");
    if (p == nullptr) {
      return "";
    }
    return *(static_cast<const char *const *>(p));
  }  // end of getTFELVersion

  std::string LibrariesManager::getUnitSystem(const std::string &l,
                                              const std::string &n) {
    const auto p = this->getSymbolAddress(l, n + "_unit_system");
    if (p == nullptr) {
      return "";
    }
    return *(static_cast<const char *const *>(p));
  }  // end of getUnitSystem

  unsigned short LibrariesManager::getAPIVersion(const std::string &l,
                                                 const std::string &n) {
    const auto s = n + "_api_version";
    const auto p = this->getSymbolAddress(l, s);
    if (p == nullptr) {
      return 0;
    }
    return *(static_cast<unsigned short *>(p));
  }  // end of getAPIVersion

  unsigned short LibrariesManager::getMaterialKnowledgeType(
      const std::string &l, const std::string &b) {
    return *(this->extract<unsigned short>(l, b + "_mfront_mkt"));
  }  // end of getMaterialKnowledgeType

  std::string LibrariesManager::getSource(const std::string &l,
                                          const std::string &n) {
    const auto p = this->getSymbolAddress(l, n + "_src");
    if (p == nullptr) {
      return "";
    }
    return *(static_cast<const char *const *>(p));
  }  // end of getSource

  std::string LibrariesManager::getAuthor(const std::string &l,
                                          const std::string &n) {
    const auto p = this->getSymbolAddress(l, n + "_author");
    if (p == nullptr) {
      return "";
    }
    return *(static_cast<const char *const *>(p));
  }  // end of getAuthor

  std::string LibrariesManager::getDate(const std::string &l,
                                        const std::string &n) {
    const auto p = this->getSymbolAddress(l, n + "_date");
    if (p == nullptr) {
      return "";
    }
    return *(static_cast<const char *const *>(p));
  }  // end of getDate

  std::string LibrariesManager::getValidator(const std::string &l,
                                             const std::string &n) {
    const auto p = this->getSymbolAddress(l, n + "_validator");
    if (p == nullptr) {
      return "";
    }
    return *(static_cast<const char *const *>(p));
  }  // end of getValidator

  std::string LibrariesManager::getBuildIdentifier(const std::string &l,
                                                   const std::string &n) {
    const auto p = this->getSymbolAddress(l, n + "_build_id");
    if (p == nullptr) {
      return "";
    }
    return *(static_cast<const char *const *>(p));
  }  // end of getBuildIdentifier

  std::string LibrariesManager::getInterface(const std::string &l,
                                             const std::string &n) {
    const auto p = this->getSymbolAddress(l, n + "_mfront_interface");
    if (p == nullptr) {
      return "";
    }
    return *(static_cast<const char *const *>(p));
  }  // end of getInterface

  unsigned short LibrariesManager::getBehaviourType(const std::string &l,
                                                    const std::string &b) {
    return *(this->extract<unsigned short>(l, b + "_BehaviourType"));
  }  // end of getBehaviourType

  unsigned short LibrariesManager::getBehaviourKinematic(const std::string &l,
                                                         const std::string &b) {
    return *(this->extract<unsigned short>(l, b + "_BehaviourKinematic"));
  }  // end of getBehaviourKinematic

  unsigned short LibrariesManager::getBehaviourSymmetry(const std::string &l,
                                                        const std::string &b) {
    return *(this->extract<unsigned short>(l, b + "_SymmetryType"));
  }  // end of getBehaviourSymmetry

  unsigned short LibrariesManager::getElasticStiffnessSymmetry(
      const std::string &l, const std::string &b) {
    return *(this->extract<unsigned short>(l, b + "_ElasticSymmetryType"));
  }  // end of getElasticStiffnessSymmetry

  bool LibrariesManager::requiresStiffnessTensor(const std::string &l,
                                                 const std::string &b,
                                                 const Hypothesis h) {
    const auto sn = "_requiresStiffnessTensor";
    const auto bv =
        *(this->extract<unsigned short>(l, b + "_" + toString(h) + sn, b + sn));
    return bv == 1 ? true : false;
  }  // end of requiresStiffnessTensor

  bool LibrariesManager::computesStoredEnergy(const std::string &l,
                                              const std::string &b,
                                              const Hypothesis h) {
    const auto sn = "_ComputesInternalEnergy";
    const auto bv =
        *(this->extract<unsigned short>(l, b + "_" + toString(h) + sn, b + sn));
    return bv == 1 ? true : false;
  }  // end of computesStoredEnergy

  bool LibrariesManager::computesDissipatedEnergy(const std::string &l,
                                                  const std::string &b,
                                                  const Hypothesis h) {
    const auto sn = "_ComputesDissipatedEnergy";
    const auto bv =
        *(this->extract<unsigned short>(l, b + "_" + toString(h) + sn, b + sn));
    return bv == 1 ? true : false;
  }  // end of computesDissipatedEnergy

  bool LibrariesManager::requiresThermalExpansionCoefficientTensor(
      const std::string &l, const std::string &b, const Hypothesis h) {
    const auto sn = "_requiresThermalExpansionCoefficientTensor";
    const auto bv =
        *(this->extract<unsigned short>(l, b + "_" + toString(h) + sn, b + sn));
    return bv == 1 ? true : false;
  }  // end of requiresThermalExpansionCoefficientTensor

  template <typename T>
  const T *LibrariesManager::extract(const std::string &l,
                                     const std::string &n) {
    const auto p = this->getSymbolAddress(l, n);
    if (p == nullptr) {
      raise("LibrariesManager::extract: could not load symbol '" + n + "'");
    }
    return static_cast<const T *const>(p);
  }  // end of extract

  template <typename T>
  const T *LibrariesManager::extract(const std::string &l,
                                     const std::string &n1,
                                     const std::string &n2) {
    const auto p = this->getSymbolAddress(l, n1, n2);
    if (p == nullptr) {
      raise("LibrariesManager::extract: could not load symbol '" + n1 +
            "' nor '" + n2 + "'");
    }
    return static_cast<const T *const>(p);
  }  // end of extract

  std::vector<std::string> LibrariesManager::getNames(const std::string &l,
                                                      const std::string &e,
                                                      const std::string &n) {
    std::vector<std::string> vars;
    const auto nb = *(this->extract<unsigned short>(l, e + "_n" + n));
    if (nb != 0u) {
      const auto res = this->extract<const char *const>(l, e + '_' + n);
      std::copy(res, res + nb, std::back_inserter(vars));
    }
    return vars;
  }  // end of getNames

  std::vector<std::string> LibrariesManager::getNames(const std::string &l,
                                                      const std::string &f,
                                                      const Hypothesis h,
                                                      const std::string &n) {
    std::vector<std::string> vars;
    const auto hn = toString(h);
    const auto nb = *(this->extract<unsigned short>(l, f + "_" + hn + "_n" + n,
                                                    f + "_n" + n));
    if (nb != 0u) {
      const auto res = this->extract<const char *const>(
          l, f + "_" + hn + '_' + n, f + '_' + n);
      std::copy(res, res + nb, std::back_inserter(vars));
    }
    return vars;
  }  // end of getNames

  std::vector<std::string> LibrariesManager::getGradientsNames(
      const std::string &l, const std::string &b, const Hypothesis h) {
    return this->getNames(l, b, h, "Gradients");
  }  // end of getGradientsNames

  std::vector<int> LibrariesManager::getGradientsTypes(const std::string &l,
                                                       const std::string &b,
                                                       const Hypothesis h) {
    std::vector<int> types;
    const auto hn = toString(h);
    const auto nb = *(this->extract<unsigned short>(
        l, b + "_" + hn + "_nGradients", b + "_nGradients"));
    const auto res = this->extract<const int>(
        l, b + "_" + hn + "_GradientsTypes", b + "_GradientsTypes");
    std::copy(res, res + nb, std::back_inserter(types));
    return types;
  }  // end of getGradientsTypes

  std::vector<std::string> LibrariesManager::getThermodynamicForcesNames(
      const std::string &l, const std::string &b, const Hypothesis h) {
    return this->getNames(l, b, h, "ThermodynamicForces");
  }  // end of getThermodynamicForcesNames

  std::vector<int> LibrariesManager::getThermodynamicForcesTypes(
      const std::string &l, const std::string &b, const Hypothesis h) {
    std::vector<int> types;
    const auto hn = toString(h);
    const auto nb = *(
        this->extract<unsigned short>(l, b + "_" + hn + "_nThermodynamicForces",
                                      b + "_nThermodynamicForces"));
    const auto res =
        this->extract<const int>(l, b + "_" + hn + "_ThermodynamicForcesTypes",
                                 b + "_ThermodynamicForcesTypes");
    std::copy(res, res + nb, std::back_inserter(types));
    return types;
  }  // end of getThermodynamicForcesTypes

  std::vector<std::pair<std::string, std::string>>
  LibrariesManager::getTangentOperatorBlocksNames(const std::string &l,
                                                  const std::string &b,
                                                  const Hypothesis h) {
    std::vector<std::pair<std::string, std::string>> blocks;
    const auto names = this->getNames(l, b, h, "TangentOperatorBlocks");
    const auto s = names.size();
    mgis::raise_if(
        s % 2 != 0,
        "LibrariesManager::getTangentOperatorBlocksNames: "
        "invalid declaration of the tangent operator blocks is invalid");
    for (decltype(names.size()) i = 0; i != s / 2; ++i) {
      blocks.push_back({names[2 * i], names[2 * i + 1]});
    }
    return blocks;
  }  // end of getTangentOperatorBlocksNames

  std::vector<std::string> LibrariesManager::getMaterialPropertiesNames(
      const std::string &l, const std::string &b, const Hypothesis h) {
    return this->getNames(l, b, h, "MaterialProperties");
  }  // end of getMaterialPropertiesNames

  std::vector<std::string> LibrariesManager::getInternalStateVariablesNames(
      const std::string &l, const std::string &b, const Hypothesis h) {
    return this->getNames(l, b, h, "InternalStateVariables");
  }  // end of getInternalStateVariablesNames

  std::vector<int> LibrariesManager::getInternalStateVariablesTypes(
      const std::string &l, const std::string &b, const Hypothesis h) {
    std::vector<int> types;
    const auto hn = toString(h);
    const auto nb = *(this->extract<unsigned short>(
        l, b + "_" + hn + "_nInternalStateVariables",
        b + "_nInternalStateVariables"));
    const auto res = this->extract<const int>(
        l, b + "_" + hn + "_InternalStateVariablesTypes",
        b + "_InternalStateVariablesTypes");
    std::copy(res, res + nb, std::back_inserter(types));
    return types;
  }  // end of getInternalStateVariablesTypes

  bool LibrariesManager::hasTemperatureBeenRemovedFromExternalStateVariables(
      const std::string &l, const std::string &b) {
    const auto s = b + "_TemperatureRemovedFromExternalStateVariables";
    if (!this->contains(l, s)) {
      // for backward compatibility
      return true;
    }
    return *(this->extract<unsigned short>(l, s)) == 1u;
  }  // end of hasTemperatureBeenRemovedFromExternalStateVariables

  std::vector<std::string> LibrariesManager::getExternalStateVariablesNames(
      const std::string &l, const std::string &b, const Hypothesis h) {
    return this->getNames(l, b, h, "ExternalStateVariables");
  }  // end of getMaterialPropertiesNames

  bool LibrariesManager::hasExternalStateVariablesTypes(const std::string &l,
                                                        const std::string &b,
                                                        const Hypothesis h) {
    const auto hn = toString(h);
    return (
        (this->contains(l, b + "_" + hn + "_ExternalStateVariablesTypes")) ||
        (this->contains(l, b + "_ExternalStateVariablesTypes")));
  }  // end of hasExternalStateVariablesTypes

  std::vector<int> LibrariesManager::getExternalStateVariablesTypes(
      const std::string &l, const std::string &b, const Hypothesis h) {
    std::vector<int> types;
    const auto hn = toString(h);
    const auto nb = *(this->extract<unsigned short>(
        l, b + "_" + hn + "_nExternalStateVariables",
        b + "_nExternalStateVariables"));
    const auto res = this->extract<const int>(
        l, b + "_" + hn + "_ExternalStateVariablesTypes",
        b + "_ExternalStateVariablesTypes");
    std::copy(res, res + nb, std::back_inserter(types));
    return types;
  }  // end of getExternalStateVariablesTypes

  bool LibrariesManager::contains(const std::string &l, const std::string &n) {
    return this->getSymbolAddress(l, n) != nullptr;
  }  // end of contains

  void *LibrariesManager::getSymbolAddress(const std::string &l,
                                           const std::string &n) {
    auto lib = this->loadLibrary(l);
#if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
    return reinterpret_cast<void *>(::GetProcAddress(lib, n.c_str()));
#else  /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)*/
    return ::dlsym(lib, n.c_str());
#endif /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */
  }    // end of getSymbolAddress

  void *LibrariesManager::getSymbolAddress(const std::string &l,
                                           const std::string &n1,
                                           const std::string &n2) {
    const auto p = this->getSymbolAddress(l, n1);
    if (p == nullptr) {
      return this->getSymbolAddress(l, n2);
    }
    return p;
  }  // end of getSymbolAddress

  LibrariesManager::libhandler LibrariesManager::loadLibrary(
      const std::string &l) {
    auto p = this->libraries.find(l);
    if (p == this->libraries.end()) {
      // this library has not been
      auto r = try_open(l);
      auto lib = r.first;
      if (lib == nullptr) {
        raise(
            "LibrariesManager::loadLibrary:"
            " library '" +
            l +
            "' could not be loaded, "
            "(" +
            getErrorMessage() + ")");
      }
      this->libraries.insert({l, lib});
      return lib;
    }
    return p->second;
  }  // end of loadLibrary

  void LibrariesManager::setParameter(const std::string &l,
                                      const std::string &b,
                                      const Hypothesis h,
                                      const std::string &p,
                                      const double v) {
    using fct_ptr_type = int (*)(const char *const, const double);
    const auto hn = toString(h);
    const auto ptr = getSymbolAddress(l, b + "_" + hn + "_setParameter",
                                      b + "_setParameter");
    if (ptr == nullptr) {
      mgis::raise(
          "LibrariesManager::setParameter: "
          "can't get the '" +
          b + "_setParameter' function (" + getErrorMessage() + ")");
    }
    const auto fct = reinterpret_cast<fct_ptr_type>(ptr);
    if (!fct(p.c_str(), v)) {
      mgis::raise(
          "LibrariesManager::setParameter: "
          "call to the '" +
          b + "_setParameter' function failed");
    }
  }  // end of setParameter

  void LibrariesManager::setParameter(const std::string &l,
                                      const std::string &b,
                                      const Hypothesis h,
                                      const std::string &p,
                                      const int v) {
    using fct_ptr_type = int (*)(const char *const, const int);
    const auto hn = toString(h);
    const auto ptr = getSymbolAddress(l, b + "_" + hn + "_setIntegerParameter",
                                      b + "_setIntegerParameter");
    if (ptr == nullptr) {
      mgis::raise(
          "LibrariesManager::setParameter: "
          "can't get the '" +
          b + "_setParameter' function (" + getErrorMessage() + ")");
    }
    const auto fct = reinterpret_cast<fct_ptr_type>(ptr);
    if (!fct(p.c_str(), v)) {
      mgis::raise(
          "LibrariesManager::setParameter: "
          "call to the '" +
          b + "_setParameter' function failed");
    }
  }  // end of setParameter

  void LibrariesManager::setParameter(const std::string &l,
                                      const std::string &b,
                                      const Hypothesis h,
                                      const std::string &p,
                                      const unsigned short v) {
    using fct_ptr_type = int (*)(const char *const, const unsigned short);
    const auto hn = toString(h);
    const auto ptr =
        getSymbolAddress(l, b + "_" + hn + "_setUnsignedShortParameter",
                         b + "_setUnsignedShortParameter");
    if (ptr == nullptr) {
      mgis::raise(
          "LibrariesManager::setParameter: "
          "can't get the '" +
          b + "_setParameter' function (" + getErrorMessage() + ")");
    }
    const auto fct = reinterpret_cast<fct_ptr_type>(ptr);
    if (!fct(p.c_str(), v)) {
      mgis::raise(
          "LibrariesManager::setParameter: "
          "call to the '" +
          b + "_setParameter' function failed");
    }
  }  // end of setParameter

  std::vector<std::string> LibrariesManager::getParametersNames(
      const std::string &l, const std::string &b, const Hypothesis h) {
    return this->getNames(l, b, h, "Parameters");
  }  // end of getMaterialPropertiesNames

  std::vector<int> LibrariesManager::getParametersTypes(const std::string &l,
                                                        const std::string &b,
                                                        const Hypothesis h) {
    std::vector<int> types;
    const auto hn = toString(h);
    const auto nb = *(this->extract<unsigned short>(
        l, b + "_" + hn + "_nParameters", b + "_nParameters"));
    const auto res = this->extract<const int>(
        l, b + "_" + hn + "_ParametersTypes", b + "_ParametersTypes");
    std::copy(res, res + nb, std::back_inserter(types));
    return types;
  }  // end of getInternalVariablesTypes

  static std::pair<std::string, std::string> buildSymbolsNames(
      const std::string &b, const std::string &h, const std::string &p) {
    return {b + "_" + h + "_" + p, b + "_" + p};
  }

  template <typename T>
  T LibrariesManager::getParameterDefaultValueImplementation(
      const std::string &l,
      const std::string &b,
      const Hypothesis h,
      const std::string &p) {
    const auto hn = toString(h);
    const auto pn = decomposeVariableName(p);
    const auto s1 =
        buildSymbolsNames(b, hn, pn.first + "_ParameterDefaultValue");
    if (this->contains(l, s1.first) || this->contains(l, s1.second)) {
      return *(this->extract<T>(l, s1.first, s1.second));
    }
    const auto s2 =
        buildSymbolsNames(b, hn, pn.second + "_ParameterDefaultValue");
    return *(this->extract<T>(l, s2.first, s2.second));
  }  // end of LibrariesManager_extractParameterDefaultValue

  double LibrariesManager::getParameterDefaultValue(
      const std::string &l,
      const std::string &b,
      const LibrariesManager::Hypothesis h,
      const std::string &p) {
    return this->getParameterDefaultValueImplementation<double>(l, b, h, p);
  }  // end of getParameterDefaultValue

  int LibrariesManager::getIntegerParameterDefaultValue(const std::string &l,
                                                        const std::string &b,
                                                        const Hypothesis h,
                                                        const std::string &p) {
    return this->getParameterDefaultValueImplementation<int>(l, b, h, p);
  }  // end of getIntegerParameterDefaultValue

  unsigned short LibrariesManager::getUnsignedShortParameterDefaultValue(
      const std::string &l,
      const std::string &b,
      const Hypothesis h,
      const std::string &p) {
    return this->getParameterDefaultValueImplementation<unsigned short>(l, b, h,
                                                                        p);
  }  // end of getUnsignedShortParameterDefaultValue

  bool LibrariesManager::hasBounds(const std::string &l,
                                   const std::string &b,
                                   const Hypothesis h,
                                   const std::string &n) {
    return this->hasLowerBound(l, b, h, n) || this->hasUpperBound(l, b, h, n);
  }  // end of hasBounds

  bool LibrariesManager::hasBoundImplementation(const std::string &l,
                                                const std::string &b,
                                                const Hypothesis h,
                                                const std::string &n,
                                                const std::string &bt) {
    const auto hn = toString(h);
    const auto vn = decomposeVariableName(n);
    const auto [n1, n2] = buildSymbolsNames(b, hn, vn.first + "_" + bt);
    const auto [n3, n4] = buildSymbolsNames(b, hn, vn.second + "_" + bt);
    return ((this->contains(l, n1)) || (this->contains(l, n2)) ||
            (this->contains(l, n3)) || (this->contains(l, n4)));
  }  // end of hasBoundImplementation

  long double LibrariesManager::getBoundImplementation(const std::string &l,
                                                       const std::string &b,
                                                       const Hypothesis h,
                                                       const std::string &n,
                                                       const std::string &bt) {
    const auto hn = toString(h);
    const auto vn = decomposeVariableName(n);
    const auto [n1, n2] = buildSymbolsNames(b, hn, vn.first + "_" + bt);
    if ((this->contains(l, n1)) || (this->contains(l, n2))) {
      return *(this->extract<long double>(l, n1, n2));
    }
    const auto [n3, n4] = buildSymbolsNames(b, hn, vn.second + "_" + bt);
    return *(this->extract<long double>(l, n3, n4));
  }  // end of getBoundImplementation

  bool LibrariesManager::hasLowerBound(const std::string &l,
                                       const std::string &b,
                                       const Hypothesis h,
                                       const std::string &n) {
    return this->hasBoundImplementation(l, b, h, n, "LowerBound");
  }  // end of hasLowerBound

  bool LibrariesManager::hasUpperBound(const std::string &l,
                                       const std::string &b,
                                       const Hypothesis h,
                                       const std::string &n) {
    return this->hasBoundImplementation(l, b, h, n, "UpperBound");
  }  // end of hasUpperBound

  long double LibrariesManager::getLowerBound(const std::string &l,
                                              const std::string &b,
                                              const Hypothesis h,
                                              const std::string &n) {
    return this->getBoundImplementation(l, b, h, n, "LowerBound");
  }  // end of getLowerBound

  long double LibrariesManager::getUpperBound(const std::string &l,
                                              const std::string &b,
                                              const Hypothesis h,
                                              const std::string &n) {
    return this->getBoundImplementation(l, b, h, n, "UpperBound");
  }  // end of getUpperBound

  bool LibrariesManager::hasPhysicalBounds(const std::string &l,
                                           const std::string &b,
                                           const Hypothesis h,
                                           const std::string &n) {
    return this->hasLowerPhysicalBound(l, b, h, n) ||
           this->hasUpperPhysicalBound(l, b, h, n);
  }  // end of hasPhysicalBounds

  bool LibrariesManager::hasLowerPhysicalBound(const std::string &l,
                                               const std::string &b,
                                               const Hypothesis h,
                                               const std::string &n) {
    return this->hasBoundImplementation(l, b, h, n, "LowerPhysicalBound");
  }  // end of hasLowerPhysicalBound

  bool LibrariesManager::hasUpperPhysicalBound(const std::string &l,
                                               const std::string &b,
                                               const Hypothesis h,
                                               const std::string &n) {
    return this->hasBoundImplementation(l, b, h, n, "UpperPhysicalBound");
  }  // end of hasUpperPhysicalBound

  long double LibrariesManager::getLowerPhysicalBound(const std::string &l,
                                                      const std::string &b,
                                                      const Hypothesis h,
                                                      const std::string &n) {
    return this->getBoundImplementation(l, b, h, n, "LowerPhysicalBound");
  }  // end of getLowerPhysicalBound

  long double LibrariesManager::getUpperPhysicalBound(const std::string &l,
                                                      const std::string &b,
                                                      const Hypothesis h,
                                                      const std::string &n) {
    return this->getBoundImplementation(l, b, h, n, "UpperPhysicalBound");
  }  // end of getUpperPhysicalBound

  LibrariesManager::~LibrariesManager() {
    for (const auto &l : this->libraries) {
#if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
      ::FreeLibrary(l.second);
#else  /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */
      ::dlclose(l.second);
#endif /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */
    }
  }  // end of ~LibrariesManager

}  // namespace mgis
