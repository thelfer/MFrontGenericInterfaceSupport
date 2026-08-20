/*!
 * \file   include/MGIS/Profiling.hxx
 * \brief
 * \author Julien Rigal, Raphaël Prat
 * \date   01/08/2026
 * \copyright (C) Copyright Thomas Helfer 2018.
 * Use, modification and distribution are subject
 * to one of the following licences:
 * - GNU Lesser General Public License (LGPL), Version 3.0. (See accompanying
 *   file LGPL-3.0.txt)
 * - CECILL-C,  Version 1.0 (See accompanying files
 *   CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt).
 */

#ifndef LIB_MGIS_PROFILING_HXX
#define LIB_MGIS_PROFILING_HXX 1

#include <chrono>
#include <string>
#include "MGIS/Config.hxx"

#define MGIS_CONCAT(a, b) MGIS_CONCAT_INNER(a, b)
#define MGIS_CONCAT_INNER(a, b) a##b
#define MGIS_VARNAME() MGIS_CONCAT(mgis_timer_, __COUNTER__)

/*Temporary macro to allow compiling (careful, the old CatchTimeSection(NAME) will not work)*/
#define MGIS_CATCH_TIME_SECTION_1(NAME) \
  mgis::ProfilingSection MGIS_VARNAME()(NAME)

#define MGIS_GET_CATCH_MACRO(_1, _2, MACRO_NAME, ...) MACRO_NAME
#define MGIS_EXPAND(x) x

#define CatchTimeSection(...) \
  MGIS_EXPAND(MGIS_GET_CATCH_MACRO(__VA_ARGS__, MGIS_CATCH_TIME_SECTION_2, MGIS_CATCH_TIME_SECTION_1)(__VA_ARGS__))

#define CatchTimeSection(CTX, NAME) \
  mgis::ProfilingSection MGIS_VARNAME()(CTX, NAME, (CTX).isProfilingEnabled())
#define CatchLocalTimeSection(CTX, NAME, IS_ENABLED) \
  mgis::ProfilingSection MGIS_VARNAME()(CTX, NAME, IS_ENABLED)



namespace mgis {

  struct Context;

  struct MGIS_EXPORT ProfilingSection {
    //! \brief Standard constructor (active or inactive depending on the 'enabled' flag)
    ProfilingSection(Context& ctx,
                     std::string name,
                     bool enabled) noexcept;
    
    //! \brief dummy constructor for 1-argument CatchTimeSection
    ProfilingSection(std::string name) noexcept;
    
    //! \brief Default constructor (fallback, always inactive)
    ProfilingSection() noexcept : ctx_ptr(nullptr), active(false) {}

    ~ProfilingSection() noexcept;

    ProfilingSection(const ProfilingSection&) = delete;
    ProfilingSection& operator=(const ProfilingSection&) = delete;

    ProfilingSection(ProfilingSection&&) = delete;
    ProfilingSection& operator=(ProfilingSection&&) = delete;

   private:
    Context* ctx_ptr;
    bool active;
    std::chrono::high_resolution_clock::time_point start;
  };

}  // end of namespace mgis

#endif
