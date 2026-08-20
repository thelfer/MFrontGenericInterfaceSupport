/*!
 * \file   include/MGIS/ProfilingData.hxx
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

}  // end of namespace mgis

#endif
