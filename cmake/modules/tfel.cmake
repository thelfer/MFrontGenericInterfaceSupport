# type of architecture
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  add_definitions("-DTFEL_ARCH64")
else(CMAKE_SIZEOF_VOID_P EQUAL 8)
  add_definitions("-DTFEL_ARCH32")
endif(CMAKE_SIZEOF_VOID_P EQUAL 8)

message(STATUS "tfel version          : ${TFEL_VERSION}")
message(STATUS "tfel C++ standard     : ${TFEL_CXX_STANDARD}")
message(STATUS "mfront                : ${MFRONT}")
message(STATUS "tfel-config           : ${TFEL_CONFIG}")
if(TFEL_CHECK)
  message(STATUS "tfel-check            : ${TFEL_CHECK}")
endif(TFEL_CHECK)
message(STATUS "tfel include          : ${TFEL_INCLUDE_PATH}")
message(STATUS "tfel libs             : ${TFEL_LIBRARIES}")
message(STATUS "tfel python bindings  : ${TFEL_PYTHON_VERSION}")
message(STATUS "TFELTests             : ${TFELTests_LIBRARY}")
message(STATUS "TFELException         : ${TFELException_LIBRARY}")
message(STATUS "TFELUtilities         : ${TFELUtilities_LIBRARY}")
message(STATUS "TFELMath              : ${TFELMath_LIBRARY}")
message(STATUS "TFELMaterial          : ${TFELMaterial_LIBRARY}")
message(STATUS "TFELPhysicalConstants : ${TFELPhysicalConstants_LIBRARY}")
message(STATUS "MTestFileGenerator    : ${MTestFileGenerator_LIBRARY}")

set(TFEL_MFRONT_LIBRARIES
    ${TFEL_LIBRARIES}
    CACHE INTERNAL "")
set(MGIS_HAVE_TFEL ON)
