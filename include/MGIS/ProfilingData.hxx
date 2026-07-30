#ifndef LIB_MGIS_PROFILING_DATA_HXX
#define LIB_MGIS_PROFILING_DATA_HXX

#include <string>
#include <vector>
#include <memory>
#include <chrono>

namespace mgis {

  struct ProfilingData {
    std::string name;
    double time_in_seconds = 0.0;
    std::size_t calls = 0;
    // For building profiling tree
    std::vector<std::unique_ptr<ProfilingData>> children;
  };

} // end of namespace mgis

#endif