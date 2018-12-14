add_library(MFrontGenericInterface-FEniCS SHARED
            QuadratureFunction.cxx)
target_include_directories(MFrontGenericInterface-FEniCS
  SYSTEM ${DOLFIN_INCLUDE_DIRS} ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})
target_link_libraries(MFrontGenericInterface-FEniCS ${DOLFIN_LIBRARIES})
