/*!
 * \file   src/VerbosityLevel.cxx
 * \brief  This file implements function related to the `VerbosityLevel`
 * enumeration \date   08/02/2023
 */

#include "MGIS/Context.hxx"
#include "MGIS/VerbosityLevel.hxx"

namespace mgis {

  VerbosityLevel &getDefaultVerbosityLevel() noexcept {
    //   if (Global::debug) {
    //     static VerbosityLevel verboseMode = verboseDebug;
    //     return verboseMode;
    //   }
    //   else {
    static VerbosityLevel verboseMode = verboseLevel1;
    return verboseMode;
    //  }
  }  // end of getDefaultVerbosityLevel()

  void setDefaultVerbosityLevel(const VerbosityLevel l) noexcept {
    getDefaultVerbosityLevel() = l;
  }  // end of setDefaultVerbosityLevel

  bool setDefaultVerbosityLevel(Context &ctx, std::string_view l) noexcept {
    const auto lvl = convertToVerbosityLevel(ctx, l);
    if (lvl.has_value()) {
      setDefaultVerbosityLevel(*lvl);
      ctx.setVerbosityLevel(*lvl);
      return true;
    }
    return false;
  }  // end of setDefaultVerbosityLevel

  std::optional<VerbosityLevel> convertToVerbosityLevel(
      Context &ctx, std::string_view l) noexcept {
    if (l == "quiet") {
      return verboseQuiet;
    } else if (l == "level0") {
      return verboseLevel0;
    } else if (l == "level1") {
      return verboseLevel1;
    } else if (l == "level2") {
      return verboseLevel2;
    } else if (l == "level3") {
      return verboseLevel3;
    } else if (l == "debug") {
      return verboseDebug;
    } else if (l == "full") {
      return verboseFull;
    }
    return ctx.registerErrorMessage(
        "convertToVerbosityLevel: unsupported verbose level '" +
        std::string{l} + "'");
  }  // end of convertToVerbosityLevel

}  // namespace mgis
