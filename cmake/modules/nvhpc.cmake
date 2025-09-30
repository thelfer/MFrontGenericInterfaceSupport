if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
  set(OPTIMISATION_FLAGS "-O2 -DNDEBUG ${OPTIMISATION_FLAGS}")
  mgis_enable_cxx_compiler_flag(OPTIMISATION_FLAGS_MARCH "fast")
endif((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions("-g")
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

add_definitions("-Ktrap=fp -Kieee -noswitcherror -Wno-long-long -Wall -Wshadow -Wextra -Wno-unused-parameter -Werror --display_error_number --diag_suppress1 --diag_suppress185 --diag_suppress941 -Wunused-variable --diag_suppress177 --diag_suppress191 --diag_suppress3312 --diag_suppress2811")

if(enable-fast-math)
  mgis_enable_cxx_compiler_flag(OPTIMISATION_FLAGS  "ffinite-math-only")
else(enable-fast-math)
  mgis_enable_cxx_compiler_flag(OPTIMISATION_FLAGS2 "ffinite-math-only")
endif(enable-fast-math)

if(enable-gpu-offloading)
  add_compile_options(-stdpar)
  list(APPEND MGIS_ADDITIONAL_LINK_FLAGS "-stdpar")
endif(enable-gpu-offloading)