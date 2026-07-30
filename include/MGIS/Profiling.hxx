#ifndef LIB_MGIS_PROFILING_HXX
#define LIB_MGIS_PROFILING_HXX 1

#include <chrono>
#include <string>
#include "MGIS/Config.hxx"

#define MGIS_CONCAT(a, b) MGIS_CONCAT_INNER(a, b)
#define MGIS_CONCAT_INNER(a, b) a##b
#define MGIS_VARNAME() MGIS_CONCAT(mgis_timer_, __COUNTER__)

#define CatchTimeSection(CTX, NAME) \
  mgis::ProfilingSection MGIS_VARNAME()(CTX, NAME, (CTX).isProfilingEnabled())
#define CatchLocalTimeSection(CTX, NAME, IS_ENABLED) \
  mgis::ProfilingSection MGIS_VARNAME()(CTX, NAME, IS_ENABLED)

namespace mgis {

  class Context;

  class MGIS_EXPORT ProfilingSection {
  public:
    //! \brief Standard constructor (active or inactive depending on the 'enabled' flag)
    ProfilingSection(Context& ctx,
                     std::string,
                     bool) noexcept;
    
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

} // end of namespace mgis

#endif