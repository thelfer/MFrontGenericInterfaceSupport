cmake_minimum_required(VERSION 3.0)

include(cmake/modules/mgis.cmake)
project("mfront-generic-interface")

# portable-build
option(enable-portable-build "produce binary that can be shared between various machine (same architecture, same gcc version, different processors" OFF)

# C-bindings (placed before compiler detection)
option(enable-c-bindings "enable c bindings support")

# fortran (placed before compiler detection)
# option(enable-fortran-bindings "enable fortran support")

# summary
#if(enable-fortran-bindings)
#  enable_language (Fortran)
#  if(NOT enable-c-bindings)
#    message(STATUS "c bindings is required for fortran bindings: "
#		   "forcing c bindings support")
#	set(enable-c-bindings ON)
#  endif(NOT enable-c-bindings)
#  message(STATUS "fortran bindings support enabled")
#endif(enable-fortran-bindings)

if(enable-c-bindings)
  message(STATUS "c bindings support enabled")
endif(enable-c-bindings)

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

#detecting TFEL and MFront
include(cmake/modules/tfel.cmake)

#documentations
add_custom_target(doc)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/docs)

# add a target to generate API documentation with Doxygen
option(enable-doxygen-doc "enable generation of the Doxygen documentation" OFF)
if(enable-doxygen-doc)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    add_custom_target(doc-html
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM)
    add_dependencies(doc doc-html)
    list(APPEND CPACK_COMPONENTS_ALL doxygen)
    set(CPACK_COMPONENT_DOXYGEN_DESCRIPTION
      "Contains the TFEL doxygen documentation")
    set(GENERATE_DOXYGEN ON)
  else(DOXYGEN_FOUND)
    message(FATAL_ERROR "doxygen is required")
  endif(DOXYGEN_FOUND)
else(enable-doxygen-doc)
  message(STATUS "doxygen documentation is not enabled")
endif(enable-doxygen-doc)

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

add_subdirectory(include)
add_subdirectory(src)
if(MGIS_HAVE_TFEL)
  add_subdirectory(tests)
endif(MGIS_HAVE_TFEL)
add_subdirectory(bindings)
