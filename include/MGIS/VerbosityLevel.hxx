/*!
 * \file   MGIS/VerbosityLevel.h
 * \brief  This file declares the `VerbosityLevel` enumeration
 * \date   08/02/2023
 */

#ifndef LIB_MGIS_VERBOSITYLEVEL_HXX
#define LIB_MGIS_VERBOSITYLEVEL_HXX 1

#include <string>
#include <optional>
#include <string_view>
#include "MGIS/Config.hxx"

namespace mgis {

  // forward declaration
  struct Context;

  /*!
   * \brief list the possible values for the logging facilities
   * provided by the `Context` class
   */
  enum VerbosityLevel {
    verboseQuiet = -1,  //<! (almost) no output
    verboseLevel0 = 0,  //<! minimal output
    verboseLevel1 = 1,  //<! a simpler output
    verboseLevel2 = 2,  //<! a much detailed output
    verboseLevel3 = 3,  //<! the finer level for standard user
    verboseDebug = 4,   //<! an output adapted for debugging
    verboseFull = 5     //<! a very detailed output
  };                    // end of enum VerboseLevel

  //! \return the current default value of the verbosity level
  [[nodiscard]] MGIS_EXPORT VerbosityLevel &getDefaultVerbosityLevel() noexcept;
  /*!
   * \brief change the default level of verbosity
   * \param[in] l: the new verbose level
   */
  MGIS_EXPORT void setDefaultVerbosityLevel(const VerbosityLevel) noexcept;
  /*!
   * \brief change the default level of verbosity
   * \param[in] ctx: execution context
   * \param[in] l: the new verbose level
   */
  MGIS_EXPORT bool setDefaultVerbosityLevel(Context &,
                                            std::string_view) noexcept;
  /*!
   * \return the level of verbosity associated with the given string
   * \param[in] ctx: execution context
   * \param[in] l: string representation of the verbose level
   */
  [[nodiscard]] MGIS_EXPORT std::optional<VerbosityLevel>
  convertToVerbosityLevel(Context &, std::string_view) noexcept;

}  // end of namespace mgis

#endif /* LIB_MGIS_VERBOSITYLEVEL_HXX */
