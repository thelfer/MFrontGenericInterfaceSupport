mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Wno-conversion")
include(cmake/modules/common-compiler-flags.cmake)

set(OPTIMISATION_FLAGS "-O2 -DNDEBUG ${OPTIMISATION_FLAGS}")
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

option(enable-sanitize-options "enable various gcc sanitize options (undefined, address,...)" OFF)

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

if(enable-sanitize-options)
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fcheck-pointer-bounds")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=bounds-strict")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=undefined")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=float-divide-by-zero")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=float-cast-overflow")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=bounds")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=alignment")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=object-size")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=vpt")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=address")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=null")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=return")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=signed-integer-overflow")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=bool")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=enum")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fstack-check")
  #  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=leak")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fno-omit-frame-pointer")
endif(enable-sanitize-options)

# basic C support
set(COMPILER_C_WARNINGS "-Wall -W -pedantic")

if(enable-parallel-stl-algorithms)
  find_package(TBB REQUIRED)
  list(APPEND MGIS_REQUIRED_ADDITIONAL_PACKAGES "TBB")
  list(APPEND MGIS_ADDITIONAL_LIBRARIES "TBB::tbb")
endif(enable-parallel-stl-algorithms)

# unsable flag
# set(COMPILER_CXXFLAGS "${COMPILER_CXXFLAGS} -D_GLIBCXX_CONCEPT_CHECKS")

if(HAVE_FORTRAN)
  include(cmake/modules/gnu-fortran-compiler.cmake)
endif(HAVE_FORTRAN)
