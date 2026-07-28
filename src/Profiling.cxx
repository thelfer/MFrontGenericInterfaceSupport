/*!
 * \file   src/Profiling.cxx
 * \brief  This file implements the `Profiling` class
 */

#include "MGIS/Profiling.hxx"
#include "MGIS/Context.hxx"

namespace mgis {

  ProfilingSection::ProfilingSection(Context& c,
                                     std::string n,
                                     bool enabled) noexcept
      : ctx_ptr(&c),
        active(enabled) {
    if (this->active && this->ctx_ptr != nullptr) {
      this->ctx_ptr->pushProfilingNode(std::move(n));
      this->start = std::chrono::high_resolution_clock::now();
    }
  } // end of ProfilingSection

  ProfilingSection::~ProfilingSection() noexcept {
    if (this->active && this->ctx_ptr != nullptr) {
      const auto end = std::chrono::high_resolution_clock::now();
      const std::chrono::duration<double> dt = end - this->start;
      
      this->ctx_ptr->popProfilingNode(dt.count());
    }
  } // end of ~ProfilingSection

} // end of namespace mgis