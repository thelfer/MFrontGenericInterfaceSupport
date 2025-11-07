if(WIN32)
  mgis_enable_cxx_compiler_flag(VISIBILITY_FLAGS "EHsc")
endif(WIN32)

mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Weverything")
mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Wno-c++98-compat-pedantic")
mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Wno-c++20-compat-pedantic")
mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Wno-padded")
mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Wno-unsafe-buffer-usage")
mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Wno-documentation")
mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Wno-documentation-unknown-command")
mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Wno-exit-time-destructors")
mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Wno-global-constructors")
mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Wno-missing-braces")
mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Wno-unsafe-buffer-usage")
mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Wrange-loop-analysis")
mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Wmove")
mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Winfinite-recursion")
mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Wcomma")
mgis_enable_cxx_compiler_flag(COMPILER_WARNINGS  "Wmicrosoft")
mgis_enable_cxx_compiler_flag2(COMPILER_WARNINGS "Wno-c++98-compat" "Wno_c__98_compat_AVAILABLE")
include(cmake/modules/common-compiler-flags.cmake)

mgis_enable_cxx_compiler_flag(OPTIMISATION_FLAGS_MARCH "march=native")
if(enable-fast-math)
  mgis_enable_cxx_compiler_flag(OPTIMISATION_FLAGS  "ffast-math")
else(enable-fast-math)
  mgis_enable_cxx_compiler_flag(OPTIMISATION_FLAGS2 "ffast-math")
endif(enable-fast-math)

if(NOT WIN32)
  mgis_enable_cxx_compiler_flag(VISIBILITY_FLAGS "fvisibility=hidden")
  mgis_enable_cxx_compiler_flag(VISIBILITY_FLAGS "fvisibility-inlines-hidden")
  set(COMPILER_DEFAULT_VISIBILITY_FLAG "-fvisibility=default")
endif(NOT WIN32)

set(OPTIMISATION_FLAGS "-DNO_RUNTIME_CHECK_BOUNDS ${OPTIMISATION_FLAGS}")

if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
  set(OPTIMISATION_FLAGS "-O2 -DNDEBUG ${OPTIMISATION_FLAGS}")
endif((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))

option(enable-glibcxx-debug "use the debug version of the C++ standard as implemented by the glib C++ library" OFF)
if(enable-glibcxx-debug)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
     set(CMAKE_CXX_FLAGS_DEBUG "-g -D_GLIBCXX_DEBUG -Rno-debug-disables-optimization" CACHE STRING
      "Flags used by the C++ compiler during debug builds."
      FORCE)
  else(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
     set(CMAKE_CXX_FLAGS_DEBUG "-g -D_GLIBCXX_DEBUG" CACHE STRING
      "Flags used by the C++ compiler during debug builds."
      FORCE)
  endif(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
else(enable-glibcxx-debug)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -Rno-debug-disables-optimization" CACHE STRING
        "Flags used by the C++ compiler during debug builds."
        FORCE)
  else(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    set(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING
        "Flags used by the C++ compiler during debug builds."
        FORCE)
  endif(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
endif(enable-glibcxx-debug)

if(HAVE_FORTRAN)
  # we associate clang with the gnu fortran compiler
  if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "LLVMFlang")
    set(LLVM_FORTRAN_COMPILER ON)
  elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM")
    set(INTEL_FORTRAN_COMPILER ON)
  else()
    include(cmake/modules/gnu-fortran-compiler.cmake)
  endif()
endif(HAVE_FORTRAN)

option(enable-libcxx "use LLVM C++ Standard library" OFF)
if(enable-libcxx)
  mgis_enable_cxx_compiler_flag(COMPILER_CXXFLAGS "stdlib=libc++")
endif(enable-libcxx)

if(enable-parallel-stl-algorithms)
 # This is a poor test to check of libstc++ is used
 if(UNIX AND NOT APPLE)
   if(NOT enable-libcxx)
     find_package(TBB REQUIRED)
     list(APPEND MGIS_REQUIRED_ADDITIONAL_PACKAGES "TBB")
     list(APPEND MGIS_ADDITIONAL_LIBRARIES "TBB::tbb")
   endif(NOT enable-libcxx)
 endif(UNIX AND NOT APPLE)
endif(enable-parallel-stl-algorithms)


option(enable-sanitize-options "enable various clang sanitize options (undefined, address,...)" OFF)

if(enable-sanitize-options)
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=address")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=thread")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=memory")
  # mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=undefined")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=cfi")
  mgis_enable_cxx_compiler_flag(COMPILER_FLAGS "fsanitize=safe-stack")
endif(enable-sanitize-options)

