mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Wno-conversion")
include(cmake/modules/common-compiler-flags.cmake)

if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release") OR
   (CMAKE_BUILD_TYPE STREQUAL "Profiling"))
  set(OPTIMISATION_FLAGS "-O2 -DNDEBUG ${OPTIMISATION_FLAGS}")
endif((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release") OR
      (CMAKE_BUILD_TYPE STREQUAL "Profiling"))
set(OPTIMISATION_FLAGS "-DNO_RUNTIME_CHECK_BOUNDS ${OPTIMISATION_FLAGS}")

mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Wdisabled-optimization")
if(NOT i586-mingw32msvc_COMPILER)
  mgis_enable_cxx_compiler_flag(VISIBILITY_FLAGS "fvisibility=hidden")
  mgis_enable_cxx_compiler_flag(VISIBILITY_FLAGS "fvisibility-inlines-hidden")
  set(COMPILER_DEFAULT_VISIBILITY_FLAG "-fvisibility=default")
endif(NOT i586-mingw32msvc_COMPILER)

mgis_enable_cxx_compiler_flag(OPTIMISATION_FLAGS_MARCH "march=native")
mgis_enable_cxx_compiler_flag(OPTIMISATION_FLAGS_MARCH "ftree-vectorize")

if (NOT CMAKE_SIZEOF_VOID_P EQUAL 8 )
  # 32 bits machines.
  # using sse and sse2 instructions rather than the
  # i387 FPU du to numerical instabilities
  mgis_enable_cxx_compiler_flag(OPTIMISATION_FLAGS_MARCH "mfpmath=sse")
  mgis_enable_cxx_compiler_flag(OPTIMISATION_FLAGS_MARCH "msse")
  mgis_enable_cxx_compiler_flag(OPTIMISATION_FLAGS_MARCH "msse2")
endif(NOT CMAKE_SIZEOF_VOID_P EQUAL 8 )

if(WIN32)
  if (CMAKE_SIZEOF_VOID_P EQUAL 8 )
    # 64 bits machines
    mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "m64")
  else(CMAKE_SIZEOF_VOID_P EQUAL 4 )
    # 32 bits machines
    mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "m32")
  endif(CMAKE_SIZEOF_VOID_P EQUAL 8 )
endif(WIN32)

if(enable-fast-math)
  mgis_enable_cxx_compiler_flag(OPTIMISATION_FLAGS  "ffast-math")
else(enable-fast-math)
  mgis_enable_cxx_compiler_flag(OPTIMISATION_FLAGS  "fno-fast-math")
  mgis_enable_cxx_compiler_flag(OPTIMISATION_FLAGS2 "ffast-math")
endif(enable-fast-math)

option(enable-glibcxx-debug "use the debug version of the C++ standard as implemented by the glib C++ library" OFF)
if(enable-glibcxx-debug)
SET(CMAKE_CXX_FLAGS_DEBUG "-g -D_GLIBCXX_DEBUG" CACHE STRING
    "Flags used by the C++ compiler during debug builds."
    FORCE)
else(enable-glibcxx-debug)
SET(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING
    "Flags used by the C++ compiler during debug builds."
    FORCE)
endif(enable-glibcxx-debug)

SET(CMAKE_C_FLAGS_DEBUG "-g" CACHE STRING
    "Flags used by the C compiler during debug builds."
    FORCE)

# coverage
SET(CMAKE_CXX_FLAGS_COVERAGE "-O0 -g -DNDEBUG -fprofile-arcs -ftest-coverage" CACHE STRING
    "Flags used by the C++ compiler during builds with tests coverage checks."
    FORCE)
SET(CMAKE_C_FLAGS_COVERAGE "-O0 -g -DNDEBUG -fprofile-arcs -ftest-coverage" CACHE STRING
    "Flags used by the C compiler during builds with tests coverage checks."
    FORCE)
set(CMAKE_EXE_LINKER_FLAGS_COVERAGE
  "-fprofile-arcs -ftest-coverage -lgcov")
set(CMAKE_MODULE_LINKER_FLAGS_COVERAGE
  "-fprofile-arcs -ftest-coverage -lgcov")
set(CMAKE_SHARED_LINKER_FLAGS_COVERAGE
  "-fprofile-arcs -ftest-coverage -lgcov")
MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_COVERAGE
  CMAKE_C_FLAGS_COVERAGE
  CMAKE_EXE_LINKER_FLAGS_COVERAGE
  CMAKE_MODULE_LINKER_FLAGS_COVERAGE
  CMAKE_SHARED_LINKER_FLAGS_COVERAGE)

# profiling
SET(CMAKE_CXX_FLAGS_PROFILING "-pg" CACHE STRING
    "Flags used by the C++ compiler during profiled builds."
    FORCE)
SET(CMAKE_C_FLAGS_PROFILING "-pg" CACHE STRING
    "Flags used by the C compiler during profiled builds."
    FORCE)
MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_PROFILING
  CMAKE_C_FLAGS_PROFILING)

# sanitizers, handled by cmake/modules/sanitizers.cmake
# note: the sanitizers available in the generic gcc documentation can be
# used: the sanitizers below are only the sanitizers supported by the
# versions of gcc tested by the CI.
# - pointer-compare and pointer-subtract require the address sanitizer
#   and the link time optimization, which are added automatically
# - hwaddress (HWASan) requires a target with memory tagging support
include(cmake/modules/sanitizers.cmake)
mgis_declare_sanitizer(address)
mgis_declare_sanitizer(thread)
mgis_declare_sanitizer(leak)
mgis_declare_sanitizer(undefined)
# sanitizers corresponding to a subset of the undefined behaviour sanitizer
# (undefined), which can be enabled individually
mgis_declare_sanitizer(alignment)
mgis_declare_sanitizer(bool)
mgis_declare_sanitizer(bounds)
mgis_declare_sanitizer(bounds-strict)
mgis_declare_sanitizer(builtin)
mgis_declare_sanitizer(enum)
mgis_declare_sanitizer(float-cast-overflow)
mgis_declare_sanitizer(float-divide-by-zero)
mgis_declare_sanitizer(null)
mgis_declare_sanitizer(object-size)
mgis_declare_sanitizer(return)
mgis_declare_sanitizer(shift)
mgis_declare_sanitizer(shift-base)
mgis_declare_sanitizer(shift-exponent)
mgis_declare_sanitizer(signed-integer-overflow)
mgis_declare_sanitizer(vptr)
# sanitizers requiring the address sanitizer and the link time optimization
mgis_declare_sanitizer(pointer-compare fsanitize=pointer-compare fsanitize=address flto)
mgis_declare_sanitizer(pointer-subtract fsanitize=pointer-subtract fsanitize=address flto)
# sanitizers requiring a target with memory tagging support
mgis_declare_sanitizer(hwaddress)
# sanitizers which are not supported by the current version of gcc: the
# dedicated options are still defined (and ignored) so that the same
# options can be used whatever the version of gcc is
mgis_declare_unavailable_sanitizer(integer)
mgis_declare_unavailable_sanitizer(vpt)
# hardening flags, which can be combined with the sanitizers above
mgis_declare_sanitizer(fstack-check fstack-check)
mgis_declare_sanitizer(fstack-protector-all fstack-protector-all)
mgis_declare_sanitizer(fno-omit-frame-pointer fno-omit-frame-pointer)
mgis_setup_sanitize_options()

# basic C support
set(COMPILER_C_WARNINGS "-Wall -W -pedantic")

if(enable-parallel-stl-algorithms)
  find_package(TBB)
  if(NOT TBB_FOUND)
     message(FATAL_ERROR "Intel Threading Building Blocks library (TBB) is required by libstdc++ to support parallel STL algorithms. You may want to disable support for those parallel algorithms by passing -Denable-parallel-stl-algorithms=OFF to cmake")
  endif(NOT TBB_FOUND)
  list(APPEND MGIS_REQUIRED_ADDITIONAL_PACKAGES "TBB")
  list(APPEND MGIS_ADDITIONAL_LIBRARIES "TBB::tbb")
endif(enable-parallel-stl-algorithms)

# unsable flag
# set(COMPILER_CXXFLAGS "${COMPILER_CXXFLAGS} -D_GLIBCXX_CONCEPT_CHECKS")

if(HAVE_FORTRAN)
  include(cmake/modules/gnu-fortran-compiler.cmake)
endif(HAVE_FORTRAN)
