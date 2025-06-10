if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
  set(OPTIMISATION_FLAGS "-O2 -DNDEBUG ${OPTIMISATION_FLAGS}")
  tfel_enable_cxx_compiler_flag(OPTIMISATION_FLAGS_MARCH "fast")
endif((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions("-g")
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

add_definitions("-fPIC -Ktrap=fp -Kieee -noswitcherror -std=c++17 -Wno-long-long -Wall -Wshadow -Wextra -Wno-unused-parameter -Werror --display_error_number --diag_suppress1 --diag_suppress185 --diag_suppress941 -Wunused-variable --diag_suppress177 --diag_suppress191")

if(enable-fast-math)
  tfel_enable_cxx_compiler_flag(OPTIMISATION_FLAGS  "ffinite-math-only")
else(enable-fast-math)
  tfel_enable_cxx_compiler_flag(OPTIMISATION_FLAGS2 "ffinite-math-only")
endif(enable-fast-math)
