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
  }    // end of  getErrorMessage

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

  static std::string decomposeVariableName(const std::string &n) {
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
  }  // end of decomposeVariableName

  LibrariesManager &LibrariesManager::get() {
    static LibrariesManager lm;
    return lm;
  }  // end of LibrariesManager::get

  LibrariesManager::LibrariesManager() = default;

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
  }  // end of LibrariesManager::getBehaviour

  std::string LibrariesManager::getTFELVersion(const std::string &l,
                                               const std::string &n) {
    const auto p = this->getSymbolAddress(l, n + "_tfel_version");
    if (p == nullptr) {
      return "";
    }
    return *(static_cast<const char *const *>(p));
  }  // end of LibrariesManager::getTFELVersion

  unsigned short LibrariesManager::getMaterialKnowledgeType(
      const std::string &l, const std::string &b) {
    return *(this->extract<unsigned short>(l, b + "_mfront_mkt"));
  }  // end of LibrariesManager::getMaterialKnowledgeType

  std::string LibrariesManager::getSource(const std::string &l,
                                          const std::string &n) {
    const auto p = this->getSymbolAddress(l, n + "_src");
    if (p == nullptr) {
      return "";
    }
    return *(static_cast<const char *const *>(p));
  }  // end of LibrariesManager::getSource

  std::string LibrariesManager::getInterface(const std::string &l,
                                             const std::string &n) {
    const auto p = this->getSymbolAddress(l, n + "_mfront_interface");
    if (p == nullptr) {
      return "";
    }
    return *(static_cast<const char *const *>(p));
  }  // end of LibrariesManager::getInterface

  unsigned short LibrariesManager::getBehaviourType(const std::string &l,
                                                    const std::string &b) {
    return *(this->extract<unsigned short>(l, b + "_BehaviourType"));
  }  // end of LibrariesManager::getBehaviourType

  unsigned short LibrariesManager::getBehaviourKinematic(const std::string &l,
                                                         const std::string &b) {
    return *(this->extract<unsigned short>(l, b + "_BehaviourKinematic"));
  }  // end of LibrariesManager::getBehaviourKinematic

  unsigned short LibrariesManager::getBehaviourSymmetry(const std::string &l,
                                                        const std::string &b) {
    return *(this->extract<unsigned short>(l, b + "_SymmetryType"));
  }  // end of LibrariesManager::getBehaviourSymmetry

  unsigned short LibrariesManager::getElasticStiffnessSymmetry(
      const std::string &l, const std::string &b) {
    return *(this->extract<unsigned short>(l, b + "_ElasticSymmetryType"));
  }  // end of LibrariesManager::getElasticStiffnessSymmetry

  bool LibrariesManager::requiresStiffnessTensor(const std::string &l,
                                                 const std::string &b,
                                                 const Hypothesis h) {
    const auto sn = "_requiresStiffnessTensor";
    const auto bv =
        *(this->extract<unsigned short>(l, b + "_" + toString(h) + sn, b + sn));
    return bv == 1 ? true : false;
  }  // end of LibrariesManager::requiresStiffnessTensor

  bool LibrariesManager::requiresThermalExpansionCoefficientTensor(
      const std::string &l, const std::string &b, const Hypothesis h) {
    const auto sn = "_requiresThermalExpansionCoefficientTensor";
    const auto bv =
        *(this->extract<unsigned short>(l, b + "_" + toString(h) + sn, b + sn));
    return bv == 1 ? true : false;
  }  // end of LibrariesManager::requiresThermalExpansionCoefficientTensor

  template <typename T>
  const T *LibrariesManager::extract(const std::string &l,
                                     const std::string &n) {
    const auto p = this->getSymbolAddress(l, n);
    if (p == nullptr) {
      raise("LibrariesManager::extract: could not load symbol '" + n + "'");
    }
    return static_cast<const T *const>(p);
  }  // end of LibrariesManager::extract

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
  }  // end of LibrariesManager::extract

  std::vector<std::string> LibrariesManager::getNames(const std::string &l,
                                                      const std::string &f,
                                                      const Hypothesis h,
                                                      const std::string &n) {
    std::vector<std::string> vars;
    const auto hn = toString(h);
    const auto nb = *(this->extract<unsigned short>(l, f + "_" + hn + "_n" + n,
                                                    f + "_n" + n));
    const auto res = this->extract<const char *const>(l, f + "_" + hn + '_' + n,
                                                      f + '_' + n);
    std::copy(res, res + nb, std::back_inserter(vars));
    return vars;
  }  // end of LibrariesManager::getNames

  std::vector<std::string> LibrariesManager::getGradientsNames(
      const std::string &l, const std::string &b, const Hypothesis h) {
    return this->getNames(l, b, h, "Gradients");
  }  // end of LibrariesManager::getGradientsNames

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
  }  // end of LibrariesManager::getGradientsTypes

  std::vector<std::string> LibrariesManager::getThermodynamicForcesNames(
      const std::string &l, const std::string &b, const Hypothesis h) {
    return this->getNames(l, b, h, "ThermodynamicForces");
  }  // end of LibrariesManager::getThermodynamicForcesNames

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
  }  // end of LibrariesManager::getThermodynamicForcesTypes

  std::vector<std::pair<std::string, std::string>>
  LibrariesManager::getTangentOperatorBlocksNames(const std::string &l,
                                                  const std::string &b,
                                                  const Hypothesis h) {
    std::vector<std::pair<std::string, std::string>> blocks;
    const auto names = this->getNames(l, b, h, "TangentOperatorBlocks");
    const auto s = names.size();
    mgis::raise_if(
        s % 2 != 0,
        "LibrariesManager::getUMATTangentOperatorBlocksNames: "
        "invalid declaration of the tangent operator blocks is invalid");
    for (decltype(names.size()) i = 0; i != s / 2; ++i) {
      blocks.push_back({names[2 * i], names[2 * i + 1]});
    }
    return blocks;
  }  // end of LibrariesManager::getTangentOperatorBlocksNames

  std::vector<std::string> LibrariesManager::getMaterialPropertiesNames(
      const std::string &l, const std::string &b, const Hypothesis h) {
    return this->getNames(l, b, h, "MaterialProperties");
  }  // end of LibrariesManager::getMaterialPropertiesNames

  std::vector<std::string> LibrariesManager::getInternalStateVariablesNames(
      const std::string &l, const std::string &b, const Hypothesis h) {
    return this->getNames(l, b, h, "InternalStateVariables");
  }  // end of LibrariesManager::getInternalStateVariablesNames

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
  }  // end of LibrariesManager::getInternalStateVariablesTypes

  std::vector<std::string> LibrariesManager::getExternalStateVariablesNames(
      const std::string &l, const std::string &b, const Hypothesis h) {
    return this->getNames(l, b, h, "ExternalStateVariables");
  }  // end of LibrariesManager::getMaterialPropertiesNames

  bool LibrariesManager::contains(const std::string &l, const std::string &n) {
    return this->getSymbolAddress(l, n) != nullptr;
  }  // end of LibrariesManager::contains

  void *LibrariesManager::getSymbolAddress(const std::string &l,
                                           const std::string &n) {
    auto lib = this->loadLibrary(l);
#if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
    return reinterpret_cast<void*>(::GetProcAddress(lib, n.c_str()));
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
  }  // end of LibrariesManager::loadLibrary

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
  }  // end of LibrariesManager::setParameter

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
  }  // end of LibrariesManager::setParameter

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
  }  // end of LibrariesManager::setParameter

  std::vector<std::string> LibrariesManager::getParametersNames(
      const std::string &l, const std::string &b, const Hypothesis h) {
    return this->getNames(l, b, h, "Parameters");
  }  // end of LibrariesManager::getMaterialPropertiesNames

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
  }  // end of LibrariesManager::getInternalVariablesTypes

  double LibrariesManager::getParameterDefaultValue(const std::string &l,
                                                    const std::string &b,
                                                    const Hypothesis h,
                                                    const std::string &p) {
    const auto pn = decomposeVariableName(p);
    const auto hn = toString(h);
    return *(this->extract<double>(
        l, b + "_" + hn + "_" + pn + "_ParameterDefaultValue",
        b + "_" + pn + "_ParameterDefaultValue"));
  }  // end of LibrariesManager::getParameterDefaultValue

  int LibrariesManager::getIntegerParameterDefaultValue(const std::string &l,
                                                        const std::string &b,
                                                        const Hypothesis h,
                                                        const std::string &p) {
    const auto pn = decomposeVariableName(p);
    const auto hn = toString(h);
    return *(this->extract<int>(
        l, b + "_" + hn + "_" + pn + "_ParameterDefaultValue",
        b + "_" + pn + "_ParameterDefaultValue"));
  }  // end of LibrariesManager::getIntegerParameterDefaultValue

  unsigned short LibrariesManager::getUnsignedShortParameterDefaultValue(
      const std::string &l,
      const std::string &b,
      const Hypothesis h,
      const std::string &p) {
    const auto pn = decomposeVariableName(p);
    const auto hn = toString(h);
    return *(this->extract<unsigned short>(
        l, b + "_" + hn + "_" + pn + "_ParameterDefaultValue",
        b + "_" + pn + "_ParameterDefaultValue"));
  }  // end of LibrariesManager::getUnsignedShortParameterDefaultValue

  bool LibrariesManager::hasBounds(const std::string &l,
                                   const std::string &b,
                                   const Hypothesis h,
                                   const std::string &n) {
    const auto hn = toString(h);
    const auto vn = decomposeVariableName(n);
    const auto n1 = b + "_" + hn + "_" + vn + "_LowerBound";
    const auto n2 = b + "_" + hn + "_" + vn + "_UpperBound";
    const auto n3 = b + "_" + vn + "_LowerBound";
    const auto n4 = b + "_" + vn + "_UpperBound";
    return ((this->contains(l, n1)) || (this->contains(l, n2)) ||
            (this->contains(l, n3)) || (this->contains(l, n4)));
  }  // end of LibrariesManager::hasBounds

  bool LibrariesManager::hasLowerBound(const std::string &l,
                                       const std::string &b,
                                       const Hypothesis h,
                                       const std::string &n) {
    const auto hn = toString(h);
    const auto vn = decomposeVariableName(n);
    const auto n1 = b + "_" + hn + "_" + vn + "_LowerBound";
    const auto n2 = b + "_" + vn + "_LowerBound";
    return ((this->contains(l, n1)) || (this->contains(l, n2)));
  }  // end of LibrariesManager::hasLowerBound

  bool LibrariesManager::hasUpperBound(const std::string &l,
                                       const std::string &b,
                                       const Hypothesis h,
                                       const std::string &n) {
    const auto hn = toString(h);
    const auto vn = decomposeVariableName(n);
    const auto n1 = b + "_" + hn + "_" + vn + "_UpperBound";
    const auto n2 = b + "_" + vn + "_UpperBound";
    return ((this->contains(l, n1)) || (this->contains(l, n2)));
  }  // end of LibrariesManager::hasUpperBound

  long double LibrariesManager::getLowerBound(const std::string &l,
                                              const std::string &b,
                                              const Hypothesis h,
                                              const std::string &n) {
    const auto hn = toString(h);
    const auto vn = decomposeVariableName(n);
    const auto n1 = b + "_" + hn + "_" + vn + "_LowerBound";
    const auto n2 = b + "_" + vn + "_LowerBound";
    return *(this->extract<long double>(l, n1, n2));
  }  // end of LibrariesManager::getLowerBound

  long double LibrariesManager::getUpperBound(const std::string &l,
                                              const std::string &b,
                                              const Hypothesis h,
                                              const std::string &n) {
    const auto hn = toString(h);
    const auto vn = decomposeVariableName(n);
    const auto n1 = b + "_" + hn + "_" + vn + "_UpperBound";
    const auto n2 = b + "_" + vn + "_UpperBound";
    return *(this->extract<long double>(l, n1, n2));
  }  // end of LibrariesManager::getUpperBound

  bool LibrariesManager::hasPhysicalBounds(const std::string &l,
                                           const std::string &b,
                                           const Hypothesis h,
                                           const std::string &n) {
    const auto hn = toString(h);
    const auto vn = decomposeVariableName(n);
    const auto n1 = b + "_" + hn + "_" + vn + "_LowerPhysicalBound";
    const auto n2 = b + "_" + hn + "_" + vn + "_UpperPhysicalBound";
    const auto n3 = b + "_" + vn + "_LowerPhysicalBound";
    const auto n4 = b + "_" + vn + "_UpperPhysicalBound";
    return ((this->contains(l, n1)) || (this->contains(l, n2)) ||
            (this->contains(l, n3)) || (this->contains(l, n4)));
  }  // end of LibrariesManager::hasPhysicalBounds

  bool LibrariesManager::hasLowerPhysicalBound(const std::string &l,
                                               const std::string &b,
                                               const Hypothesis h,
                                               const std::string &n) {
    const auto hn = toString(h);
    const auto vn = decomposeVariableName(n);
    const auto n1 = b + "_" + hn + "_" + vn + "_LowerPhysicalBound";
    const auto n2 = b + "_" + vn + "_LowerPhysicalBound";
    return ((this->contains(l, n1)) || (this->contains(l, n2)));
  }  // end of LibrariesManager::hasLowerPhysicalBound

  bool LibrariesManager::hasUpperPhysicalBound(const std::string &l,
                                               const std::string &b,
                                               const Hypothesis h,
                                               const std::string &n) {
    const auto hn = toString(h);
    const auto vn = decomposeVariableName(n);
    const auto n1 = b + "_" + hn + "_" + vn + "_UpperPhysicalBound";
    const auto n2 = b + "_" + vn + "_UpperPhysicalBound";
    return ((this->contains(l, n1)) || (this->contains(l, n2)));
  }  // end of LibrariesManager::hasUpperPhysicalBound

  long double LibrariesManager::getLowerPhysicalBound(const std::string &l,
                                                      const std::string &b,
                                                      const Hypothesis h,
                                                      const std::string &n) {
    const auto hn = toString(h);
    const auto vn = decomposeVariableName(n);
    const auto n1 = b + "_" + hn + "_" + vn + "_LowerPhysicalBound";
    const auto n2 = b + "_" + vn + "_LowerPhysicalBound";
    return *(this->extract<long double>(l, n1, n2));
  }  // end of LibrariesManager::getLowerPhysicalBound

  long double LibrariesManager::getUpperPhysicalBound(const std::string &l,
                                                      const std::string &b,
                                                      const Hypothesis h,
                                                      const std::string &n) {
    const auto hn = toString(h);
    const auto vn = decomposeVariableName(n);
    const auto n1 = b + "_" + hn + "_" + vn + "_UpperPhysicalBound";
    const auto n2 = b + "_" + vn + "_UpperPhysicalBound";
    return *(this->extract<long double>(l, n1, n2));
  }  // end of LibrariesManager::getUpperPhysicalBound

  LibrariesManager::~LibrariesManager() {
    for (const auto &l : this->libraries) {
#if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
      ::FreeLibrary(l.second);
#else  /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */
      ::dlclose(l.second);
#endif /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */
    }
  }  // end of LibrariesManager::~LibrariesManager

}  // namespace mgis
