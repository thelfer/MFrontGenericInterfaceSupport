cmake_minimum_required(VERSION 3.10)
string(CONCAT MFrontGenericInterfaceSupportDescription
	"The MFrontGenericInterfaceSupport project provides "
	"useful feature for solvers' developpers to "
    "integrate MFront generated material libraries.")

include(cmake/modules/mgis.cmake)
set(MFrontGenericInterfaceSupportDevelopmentVersion ON)
project("mfront-generic-interface"
        VERSION 3.1.0
        DESCRIPTION  "${MFrontGenericInterfaceSupportDescription}"
        HOMEPAGE_URL "https://thelfer.github.io/mgis/web/index.html")

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/mgis-config-version.cmake"
  VERSION "${VERSION}"
  COMPATIBILITY SameMinorVersion)
install(FILES "${PROJECT_BINARY_DIR}/mgis-config-version.cmake"
  DESTINATION "${mgis_export_install_path}")

set(CMAKE_CXX_STANDARD 20)
set(CXX_STANDARD_REQUIRED ON)

# portable-build
option(enable-portable-build "produce binary that can be shared between various machine (same architecture, same gcc version, different processors" OFF)

# option(enable-eigen-bindings "enable eigen bindings support" OFF)
# if(enable-eigen-bindings)
#   find_package (Eigen3 REQUIRED)
# endif (enable-eigen-bindings)

# enable mgis-function
option(enable-mgis-function "enable mgis function" ON)
# enable exceptions usage
option(enable-exceptions "use exceptions to report contract violation and error reporting" OFF)

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

if(enable-python-bindings)
  find_package(Python REQUIRED COMPONENTS Interpreter Development NumPy)
  set(HAVE_PYTHON ON)
  get_filename_component(PYTHON_LIBRARY_PATH ${Python_LIBRARIES} PATH)
  get_filename_component(PYTHON_LIBRARY_FULL ${Python_LIBRARIES} NAME)
  string(REGEX REPLACE "lib(.+)\\.(.+)$" "\\1" PYTHON_LIBRARY ${PYTHON_LIBRARY_FULL})
  message(STATUS "found python          ${Python_VERSION}")
  message(STATUS "python include path   ${Python_INCLUDE_DIRS}")
  message(STATUS "python libraries path ${PYTHON_LIBRARY_PATH}")
  message(STATUS "python library        ${PYTHON_LIBRARY}")
  message(STATUS "python interpreter    ${Python_EXECUTABLE}")
  message(STATUS "Numpy version         ${Python_NumPy_VERSION}")
  message(STATUS "Numpy header path     ${Python_NumPy_INCLUDE_DIRS}")
  # pybind11
  find_package(pybind11 CONFIG REQUIRED)
  list(APPEND CPACK_COMPONENTS_ALL python_bindings)
  set(CPACK_COMPONENT_PYTHON_BINDINGS_DESCRIPTION
    "Contains the python bindings for MGIS")
  set(CPACK_COMPONENT_PYTHON_BINDINGS_DEPENDS mgis)
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
  # The current version (0.16.0) of CxxWrap overrides the
  # CMAKE_CXX_STANDARD variables, so we reset it after
  # calling find_package to the current version
  set(_CMAKE_CXX_STANDARD "${CMAKE_CXX_STANDARD}")
  find_package(JlCxx REQUIRED)
  set(CMAKE_CXX_STANDARD "${_CMAKE_CXX_STANDARD}")
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
option(disable-tfel "explicitly disable TFEL support. By default, MGIS tries to use TFEL for tests and MGIS/Function" OFF)
if(NOT disable-tfel)
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
  find_package(TFEL)
  if(TFEL_FOUND)
    include(tfel)
    include(material-properties)
    include(behaviours)
  endif(TFEL_FOUND)
endif(NOT disable-tfel)

#documentations
option(enable-doxygen-doc "enable generation of the Doxygen documentation" OFF)
option(enable-website "enable generation of the website" ON)

if(enable-doxygen-doc OR enable-website)
  add_custom_target(doc)
endif(enable-doxygen-doc OR enable-website)

# add a target to generate API documentation with Doxygen
if(enable-doxygen-doc)
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/docs)
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
      "Contains MGIS' doxygen documentation")
    set(GENERATE_DOXYGEN ON)
  else(DOXYGEN_FOUND)
    message(FATAL_ERROR "doxygen is required")
  endif(DOXYGEN_FOUND)
else(enable-doxygen-doc)
  message(STATUS "doxygen documentation is not enabled")
endif(enable-doxygen-doc)

if(enable-website)
  set(MGIS_GENERATE_WEBSITE OFF)
  # Looking for pandoc (http://johnmacfarlane.net/pandoc)
  include(cmake/modules/pandoc.cmake)
  if(MGIS_HAVE_PANDOC)
    set(MGIS_GENERATE_WEBSITE ON)
  endif(MGIS_HAVE_PANDOC)
  # Looking for jupyter-nbconvert
  include(cmake/modules/jupyter-nbconvert.cmake)
  if(MGIS_HAVE_JUPYTER_NBCONVERT)
    set(MGIS_GENERATE_WEBSITE ON)
  endif(MGIS_HAVE_JUPYTER_NBCONVERT)
  if(MGIS_GENERATE_WEBSITE)
    list(APPEND CPACK_COMPONENTS_ALL website)
    set(CPACK_COMPONENT_WEBSITE_DESCRIPTION
      "Contains the MGIS web site")
  endif(MGIS_GENERATE_WEBSITE)
else(enable-website)
  message(STATUS "Generation of the MGIS website disabled")
endif(enable-website)

add_subdirectory(docs)

# testing
set(CTEST_CONFIGURATION_TYPE "${JOB_BUILD_CONFIGURATION}")
# (must be placed *before* any add_subdirectory, cmake bug ?)
enable_testing()

if(CMAKE_VERSION VERSION_GREATER 3.21.4)
  set(ctest_args )
  list(APPEND ctest_args --output-junit ctest-results.xml) 
endif(CMAKE_VERSION VERSION_GREATER 3.21.4)
if(CMAKE_CONFIGURATION_TYPES)
  add_custom_target(check COMMAND
    ${CMAKE_CTEST_COMMAND} ${ctest_args} -T test -C $<CONFIGURATION>)
else(CMAKE_CONFIGURATION_TYPES)
  add_custom_target(check COMMAND
    ${CMAKE_CTEST_COMMAND} ${ctest_args} -T test )
endif(CMAKE_CONFIGURATION_TYPES)

add_subdirectory(include)
add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(bindings)

# Packaging
if(MFrontGenericInterfaceSupportDevelopmentVersion)
  set(CPACK_PACKAGE_VERSION_PATCH "${CMAKE_PROJECT_VERSION_PATCH}-dev")
endif()
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LGPL-3.0.txt")
set(CPACK_VERBATIM_VARIABLES YES)
if(NOT CPACK_SOURCE_GENERATOR)
  set(CPACK_SOURCE_GENERATOR "TGZ")
endif(NOT CPACK_SOURCE_GENERATOR)
set(CPACK_SOURCE_IGNORE_FILES "build;/\\.git/;.*~")
include(CPack)
