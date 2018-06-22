cmake_minimum_required(VERSION 3.0)

project("mfront-generic-interface")

# portable-build
option(enable-portable-build "produce binary that can be shared between various machine (same architecture, same gcc version, different processors" OFF)

# fortran (placed before compiler detection)
option(enable-fortran "enable fortran support")
if(enable-fortran)
  set(HAVE_FORTRAN ON)
  add_definitions("-DHAVE_FORTRAN=1")
  enable_language (Fortran)
  message(STATUS "fortran support enabled")
  set(enable-fortran-bindings ON)
else(enable-fortran)
  set(HAVE_FORTRAN OFF)
  message(STATUS "fortran support disabled")
  set(enable-fortran-bindings OFF)
endif(enable-fortran)

# Looking for libdl...
if(UNIX)
  if(NOT MFRONT_DL_LIBRARY)
    if((${CMAKE_SYSTEM_NAME} MATCHES "Linux") OR
	(${CMAKE_SYSTEM_NAME} MATCHES "Darwin"))
      find_library(MFRONT_DL_LIBRARY
	NAMES dl)
      if(NOT MFRONT_DL_LIBRARY)
	MESSAGE(FATAL_ERROR "the libdl library not found. You may want to set the `MFRONT_DL_LIBRARY` variable manually")
      endif(NOT MFRONT_DL_LIBRARY)
    else((${CMAKE_SYSTEM_NAME} MATCHES "Linux") OR
	(${CMAKE_SYSTEM_NAME} MATCHES "Darwin"))
      set(MFRONT_DL_LIBRARY "")
    endif((${CMAKE_SYSTEM_NAME} MATCHES "Linux") OR
      (${CMAKE_SYSTEM_NAME} MATCHES "Darwin"))
  endif(NOT MFRONT_DL_LIBRARY)
else(UNIX)
  set(MFRONT_DL_LIBRARY "")
endif(UNIX)

#compiler options
include(cmake/modules/compiler.cmake)
if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
  include(cmake/modules/CodeCoverage.cmake)
  SETUP_TARGET_FOR_COVERAGE(coverage ctest coverage)
endif(CMAKE_BUILD_TYPE STREQUAL "Coverage")

# testing
set(CTEST_CONFIGURATION_TYPE "${JOB_BUILD_CONFIGURATION}")
# (must be placed *before* any add_subdirectory, cmake bug ?)
enable_testing()
if(CMAKE_CONFIGURATION_TYPES)
  add_custom_target(check COMMAND 
    ${CMAKE_CTEST_COMMAND} -T test -C $<CONFIGURATION>)
else(CMAKE_CONFIGURATION_TYPES)
  add_custom_target(check COMMAND 
    ${CMAKE_CTEST_COMMAND} -T test )
endif(CMAKE_CONFIGURATION_TYPES)

add_subdirectory(src)
add_subdirectory(tests)
