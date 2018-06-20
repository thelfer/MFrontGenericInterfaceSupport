/*!
 * \file   LibrariesManager.cxx
 * \brief
 * \author Thomas Helfer
 * \date   20/06/2018
 */

#include <algorithm>
#include <cstring>
#if !((defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__))
#include <dlfcn.h>
#endif /* !((defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)) */

#include "MFront/LibrariesManager.hxx"
#include "MFront/Raise.hxx"

namespace mfront {

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
} // end of  getErrorMessage

static LibrariesManager::libhandler load_library(const std::string &l) {
#if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
  return ::LoadLibrary(TEXT(l.c_str()));
#else
  return ::dlopen(l.c_str(), RTLD_NOW);
#endif /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */
} // end of load_library

static std::pair<LibrariesManager::libhandler, std::string>
try_open(const std::string &l) {
  auto starts_with = [](const std::string &s1, const char *const s2) {
    const auto ls2 = std::strlen(s2);
    return ((s1.size() >= ls2) && (std::equal(s2, s2 + ls2, s1.begin())));
  }; // end of starts_with
  auto ends_with = [](const std::string &s1, const char *const s2) {
    const auto ls2 = std::strlen(s2);
    if (!(s1.size() >= ls2)) {
      return false;
    }
    return std::equal(s2, s2 + ls2, s1.begin() + (s1.size() - ls2));
  }; // end of ends_with
#if (defined(macintosh) || defined(Macintosh) ||                               \
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
} // end of try_open

LibrariesManager &LibrariesManager::get() {
  static LibrariesManager lm;
  return lm;
} // end of LibrariesManager::get

LibrariesManager::LibrariesManager() = default;

std::string LibrariesManager::getTFELVersion(const std::string &l,
                                             const std::string &n) {
  const auto p = this->getSymbolAddress(l, n + "_tfel_version");
  if (p == nullptr) {
    return "";
  }
  return *(static_cast<const char *const *>(p));
} // end of LibrariesManager::getTFELVersion

std::string LibrariesManager::getSource(const std::string &l,
                                        const std::string &n) {
  const auto p = this->getSymbolAddress(l, n + "_src");
  if (p == nullptr) {
    return "";
  }
  return *(static_cast<const char *const *>(p));
} // end of LibrariesManager::getSource

unsigned short LibrariesManager::getBehaviourType(const std::string &l,
                                                  const std::string &b) {
  return *(this->extract<unsigned short>(l, b + "_BehaviourType"));
} // end of LibrariesManager::getBehaviourType

unsigned short LibrariesManager::getBehaviourKinematic(const std::string &l,
                                                       const std::string &b) {
  return *(this->extract<unsigned short>(l, b + "_BehaviourKinematic"));
} // end of LibrariesManager::getBehaviourKinematic

unsigned short LibrariesManager::getBehaviourSymmetry(const std::string & l,
                                                      const std::string & b){
  return *(this->extract<unsigned short>(l, b + "_SymmetryType"));
} // end of LibrariesManager::getBehaviourSymmetry

unsigned short
LibrariesManager::getElasticStiffnessSymmetry(const std::string &l,
                                              const std::string &b) {
  return *(this->extract<unsigned short>(l, b + "_ElasticSymmetryType"));
} // end of LibrariesManager::getElasticStiffnessSymmetry

bool LibrariesManager::requiresStiffnessTensor(const std::string &l,
                                               const std::string &b,
                                               const Hypothesis h) {
  const auto sn = "_requiresStiffnessTensor";
  const auto bv =
      *(this->extract<unsigned short>(l, b + "_" + toString(h) + sn, b + sn));
  return bv == 1 ? true : false;
} // end of LibrariesManager::requiresStiffnessTensor

bool LibrariesManager::requiresThermalExpansionCoefficientTensor(
    const std::string &l, const std::string &b, const Hypothesis h) {
  const auto sn = "_requiresThermalExpansionCoefficientTensor";
  const auto bv =
      *(this->extract<unsigned short>(l, b + "_" + toString(h) + sn, b + sn));
  return bv == 1 ? true : false;
} // end of LibrariesManager::requiresThermalExpansionCoefficientTensor

template <typename T>
const T *LibrariesManager::extract(const std::string &l, const std::string &n) {
  const auto p = this->getSymbolAddress(l, n);
  if (p == nullptr) {
    raise("LibrariesManager::extract: could not load symbols '" + n + "'");
  }
  return static_cast<const T *const>(p);
} // end of LibrariesManager::extract

template <typename T>
const T *LibrariesManager::extract(const std::string &l, const std::string &n1,
                                   const std::string &n2) {
  const auto p = this->getSymbolAddress(l, n1, n2);
  if (p == nullptr) {
    raise("LibrariesManager::extract: could not load symbols '" + n1 +
          "' or '" + n1 + "'");
  }
  return static_cast<const T *const>(p);
} // end of LibrariesManager::extract

std::vector<std::string> LibrariesManager::getNames(const std::string &l,
                                                    const std::string &f,
                                                    const Hypothesis h,
                                                    const std::string &n) {
  std::vector<std::string> vars;
  const auto hn = toString(h);
  const auto nb = *(
      this->extract<unsigned short>(l, f + "_" + hn + "_n" + n, f + "_n" + n));
  const auto res =
      this->extract<const char *const>(l, f + "_" + hn + '_' + n, f + '_' + n);
  std::copy(res, res + nb, std::back_inserter(vars));
  return vars;
} // end of LibrariesManager::getNames

std::vector<std::string> LibrariesManager::getMaterialPropertiesNames(
    const std::string &l, const std::string &b, const Hypothesis h) {
  return this->getNames(l, b, h, "MaterialProperties");
} // end of LibrariesManager::getMaterialPropertiesNames

std::vector<std::string> LibrariesManager::getInternalStateVariablesNames(
    const std::string &l, const std::string &b, const Hypothesis h) {
  return this->getNames(l, b, h, "InternalStateVariables");
} // end of LibrariesManager::getMaterialPropertiesNames

std::vector<int> LibrariesManager::getInternalStateVariablesTypes(
    const std::string &l, const std::string &b, const Hypothesis h) {
  std::vector<int> types;
  const auto hn = toString(h);
  const auto nb = *(this->extract<unsigned short>(
      l, b + "_" + hn + "_nInternalStateVariables",
      b + "_nInternalStateVariables"));
  const auto res =
      this->extract<const int>(l, b + "_" + hn + "_InternalStateVariablesTypes",
                               b + "_InternalStateVariablesTypes");
  std::copy(res, res + nb, std::back_inserter(types));
  return types;
} // end of LibrariesManager::getMaterialPropertiesNames

std::vector<std::string> LibrariesManager::getExternalStateVariablesNames(
    const std::string &l, const std::string &b, const Hypothesis h) {
  return this->getNames(l, b, h, "ExternalStateVariables");
} // end of LibrariesManager::getMaterialPropertiesNames


void *LibrariesManager::getSymbolAddress(const std::string &l,
                                         const std::string &n) {
  auto lib = this->loadLibrary(l);
#if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
  return ::GetProcAddress(lib, n.c_str());
#else  /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)*/
  return ::dlsym(lib, n.c_str());
#endif /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */
} // end of getSymbolAddress

void *LibrariesManager::getSymbolAddress(const std::string &l,
                                         const std::string &n1,
                                         const std::string &n2) {
  const auto p = this->getSymbolAddress(l, n1);
  if (p == nullptr) {
    return this->getSymbolAddress(l, n2);
  }
  return p;
} // end of getSymbolAddress

LibrariesManager::libhandler
LibrariesManager::loadLibrary(const std::string &l) {
  auto p = this->libraries.find(l);
  if (p == this->libraries.end()) {
    // this library has not been
    auto r = try_open(l);
    auto lib = r.first;
    if (lib == nullptr) {
      raise("LibrariesManager::loadLibrary:"
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
} // end of LibrariesManager::loadLibrary

LibrariesManager::~LibrariesManager() {
  for (const auto &l : this->libraries) {
#if (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__)
    ::FreeLibrary(l.second);
#else  /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */
    ::dlclose(l.second);
#endif /* (defined _WIN32 || defined _WIN64) && (!defined __CYGWIN__) */
  }
} // end of LibrariesManager::~LibrariesManager

} // namespace mfront
