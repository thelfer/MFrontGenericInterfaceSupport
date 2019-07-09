cmake_minimum_required(VERSION 3.0)

include(cmake/modules/mgis.cmake)
project("mfront-generic-interface")

# portable-build
option(enable-portable-build "produce binary that can be shared between various machine (same architecture, same gcc version, different processors" OFF)

# C-bindings (placed before compiler detection)
option(enable-c-bindings "enable c bindings support")

# fortran-bindings (placed before compiler detection)
option(enable-fortran-bindings "enable fortran bindings support" OFF)
if(enable-fortran-bindings)
  if(NOT enable-c-bindings)
    message(FATAL_ERROR "c-bindings must be enabled for fortran bindings.")
  endif(NOT enable-c-bindings)
  enable_language (Fortran)
endif(enable-fortran-bindings)

#python bindings
option(enable-python-bindings "enable  python bindings support" OFF)
option(enable-broken-boost-python-module-visibility-handling
  "workaround a bug in boost python which occurs on FreeBSD+clang. See https://github.com/boostorg/python/issues/173 for details" OFF)

if(enable-python-bindings)
  find_package(PythonLibs REQUIRED)
  set(HAVE_PYTHON ON)
  get_filename_component(PYTHON_LIBRARY_PATH ${PYTHON_LIBRARIES} PATH)
  get_filename_component(PYTHON_LIBRARY_FULL ${PYTHON_LIBRARIES} NAME)
  string(REGEX REPLACE "lib(.+)\\.(.+)$" "\\1" PYTHON_LIBRARY ${PYTHON_LIBRARY_FULL})
  message(STATUS "found python          ${PYTHONLIBS_VERSION_STRING}")
  message(STATUS "python include path   ${PYTHON_INCLUDE_DIRS}")
  message(STATUS "python libraries path ${PYTHON_LIBRARY_PATH}")
  message(STATUS "python library        ${PYTHON_LIBRARY}")
  string(REGEX REPLACE "[a-z]+.*$" "" PYTHONLIBS_VERSION_CLEANED ${PYTHONLIBS_VERSION_STRING})
  find_package(PythonInterp ${PYTHONLIBS_VERSION_CLEANED} REQUIRED)
  message(STATUS "python interpreter    ${PYTHON_EXECUTABLE}")
  if(NOT HAVE_PYTHON)
    message(FATAL_ERROR "python not found. python is mandatory for python bindings.")
  endif(NOT HAVE_PYTHON)
  # numpy detection
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c
            "from __future__ import print_function\ntry: import numpy; print(numpy.get_include(), end='')\nexcept:pass\n"
            OUTPUT_VARIABLE NUMPY_INCLUDE_PATH)
  find_path(PYTHON_NUMPY_INCLUDE_DIR numpy/arrayobject.h
    HINTS "${NUMPY_INCLUDE_PATH}" "${PYTHON_INCLUDE_PATH}" NO_DEFAULT_PATH)
  if(NOT PYTHON_NUMPY_INCLUDE_DIR)
    message(FATAL_ERROR "numpy not found. numpy is mandatory for python bindings.")
  endif(NOT PYTHON_NUMPY_INCLUDE_DIR)
  # boost python detection
  message(STATUS "Trying to find libboost_python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
  find_package(Boost 1.36.0 COMPONENTS
               "python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
  if(Boost_PYTHON-PY${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}_FOUND)
    set(Boost_PYTHON_FOUND ON)
    set(Boost_PYTHON_LIBRARY "${Boost_PYTHON-PY${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}_LIBRARY}")
  else(Boost_PYTHON-PY${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}_FOUND)
    set(Boost_PYTHON_FOUND OFF)
  endif(Boost_PYTHON-PY${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}_FOUND)
  if(NOT Boost_PYTHON_FOUND)
	  message(STATUS "Trying to find libboost_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
	  find_package(Boost 1.36.0 COMPONENTS
                   "python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
	  if(Boost_PYTHON${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}_FOUND)
	    set(Boost_PYTHON_FOUND ON)
	    set(Boost_PYTHON_LIBRARY "${Boost_PYTHON${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}_LIBRARY}")
	  else(Boost_PYTHON${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}_FOUND)
	    set(Boost_PYTHON_FOUND OFF)
	  endif(Boost_PYTHON${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}_FOUND)
  endif(NOT Boost_PYTHON_FOUND)
  if(NOT Boost_PYTHON_FOUND)
    message(STATUS "Trying to find libboost_python-py${PYTHON_VERSION_MAJOR}")
    find_package(Boost 1.36.0 COMPONENTS
			   "python-py${PYTHON_VERSION_MAJOR}")
    if(Boost_PYTHON-PY${PYTHON_VERSION_MAJOR}_FOUND)
      set(Boost_PYTHON_FOUND ON)
      set(Boost_PYTHON_LIBRARY "${Boost_PYTHON-PY${PYTHON_VERSION_MAJOR}_LIBRARY}")
    else(Boost_PYTHON-PY${PYTHON_VERSION_MAJOR}_FOUND)
      set(Boost_PYTHON_FOUND OFF)
    endif(Boost_PYTHON-PY${PYTHON_VERSION_MAJOR}_FOUND)
  endif(NOT Boost_PYTHON_FOUND)
  if(NOT Boost_PYTHON_FOUND)
    message(STATUS "Trying to find libboost_python${PYTHON_VERSION_MAJOR}")
    find_package(Boost 1.36.0 COMPONENTS
			   "python${PYTHON_VERSION_MAJOR}")
    if(Boost_PYTHON${PYTHON_VERSION_MAJOR}_FOUND)
      set(Boost_PYTHON_FOUND ON)
      set(Boost_PYTHON_LIBRARY "${Boost_PYTHON${PYTHON_VERSION_MAJOR}_LIBRARY}")
    else(Boost_PYTHON${PYTHON_VERSION_MAJOR}_FOUND)
      set(Boost_PYTHON_FOUND OFF)
    endif(Boost_PYTHON${PYTHON_VERSION_MAJOR}_FOUND)
  endif(NOT Boost_PYTHON_FOUND)
  if(NOT Boost_PYTHON_FOUND)
    message(STATUS "Trying to find libboost_python")
    find_package(Boost 1.36.0 COMPONENTS python)
  endif(NOT Boost_PYTHON_FOUND)
  if(Boost_PYTHON_FOUND)
    message(STATUS "Boost library python found (${Boost_PYTHON_LIBRARY}).")
  else(Boost_PYTHON_FOUND)
    message(FATAL_ERROR "Boost python libraries not found.\n"
      "Required packages are :" ${BOOST_MANDATORY_COMPONENTS})
  endif(Boost_PYTHON_FOUND)
  list(APPEND CPACK_COMPONENTS_ALL python_bindings)
  set(CPACK_COMPONENT_PYTHON_BINDINGS_DESCRIPTION 
    "Contains the python bindings for TFEL, MFront and MTest")
  set(CPACK_COMPONENT_PYTHON_BINDINGS_DEPENDS  core mfront mtest)
  if(enable-broken-boost-python-module-visibility-handling)
      add_definitions("-DMGIS_BOOST_PYTHON_USE_GCC_SYMBOL_VISIBILITY")
  endif(enable-broken-boost-python-module-visibility-handling)
endif(enable-python-bindings)

# fenics-bindings (placed before compiler detection)
option(enable-fenics-bindings "enable fenics bindings support" OFF)
if(enable-fenics-bindings)
  # Find DOLFIN
  find_package(DOLFIN REQUIRED)
  include(${DOLFIN_USE_FILE})
endif(enable-fenics-bindings)

# julia-bindings (placed before compiler detection)
option(enable-julia-bindings "enable julia bindings support" OFF)
if(enable-julia-bindings)
  find_package(JlCxx REQUIRED)
endif(enable-julia-bindings)


# summary

if(enable-c-bindings)
  message(STATUS "c bindings support enabled")
endif(enable-c-bindings)

if(enable-python-bindings)
  message(STATUS "python bindings support enabled")
endif(enable-python-bindings)

if(enable-fortran-bindings)
  message(STATUS "fortran bindings support enabled")
endif(enable-fortran-bindings)

if(enable-fenics-bindings)
  message(STATUS "FEniCS bindings support enabled")
endif(enable-fenics-bindings)

find_package(Threads)
# if(MINGW)
#   file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/thread-mingw.cxx"
#     "int main(){std::thread t; return 0;}")
#   try_compile(MGIS_MINGW_HAVE_NATIVE_STD_THREAD
#     ${CMAKE_BINARY_DIR}
#     ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/thread-mingw.cxx
#     COMPILE_DEFINITIONS ${COMPILER_CXXFLAGS})
#   if(MGIS_MINGW_HAVE_NATIVE_STD_THREAD)
#     set(MGIS_USE_MINGW_STD_THREADS ON)
#   endif(MGIS_MINGW_HAVE_NATIVE_STD_THREAD)
# endif(MINGW)

# Looking for libdl...
if(UNIX)
  if(NOT MGIS_DL_LIBRARY)
    if((${CMAKE_SYSTEM_NAME} MATCHES "Linux") OR
	(${CMAKE_SYSTEM_NAME} MATCHES "Darwin"))
      find_library(MGIS_DL_LIBRARY
	NAMES dl)
      if(NOT MGIS_DL_LIBRARY)
	MESSAGE(FATAL_ERROR "the libdl library not found. You may want to set the `MGIS_DL_LIBRARY` variable manually")
      endif(NOT MGIS_DL_LIBRARY)
    else((${CMAKE_SYSTEM_NAME} MATCHES "Linux") OR
	(${CMAKE_SYSTEM_NAME} MATCHES "Darwin"))
      set(MGIS_DL_LIBRARY "")
    endif((${CMAKE_SYSTEM_NAME} MATCHES "Linux") OR
      (${CMAKE_SYSTEM_NAME} MATCHES "Darwin"))
  endif(NOT MGIS_DL_LIBRARY)
else(UNIX)
  set(MGIS_DL_LIBRARY "")
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

option(enable-website "enable generation of the website" ON)
if(enable-website)
  # Looking for pandoc (http://johnmacfarlane.net/pandoc)
  include(cmake/modules/pandoc.cmake)
  if(HAVE_PANDOC)
    list(APPEND CPACK_COMPONENTS_ALL website)
    set(CPACK_COMPONENT_WEBSITE_DESCRIPTION
      "Contains the TFEL web site")
    set(MGIS_GENERATE_WEBSITE ON)
  endif(HAVE_PANDOC)
else(enable-website)
  message(STATUS "Generation of the TFEL website disabled")
endif(enable-website)

add_subdirectory(docs)

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
