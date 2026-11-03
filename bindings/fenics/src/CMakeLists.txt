add_library(MFrontGenericInterface-FEniCS SHARED
  Utils.cxx
  NonLinearMaterial.cxx
  NonLinearMaterialFunctionBase.cxx
  NonLinearMaterialThermodynamicForcesFunction.cxx
  NonLinearMaterialTangentOperatorFunction.cxx
  NonLinearMechanicalProblem.cxx)
target_include_directories(MFrontGenericInterface-FEniCS
  PUBLIC "${PROJECT_SOURCE_DIR}/include"
  PUBLIC "${PROJECT_SOURCE_DIR}/bindings/fenics/include")
target_include_directories(MFrontGenericInterface-FEniCS
  SYSTEM PUBLIC ${DOLFIN_INCLUDE_DIRS} ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})
target_link_libraries(MFrontGenericInterface-FEniCS
  PUBLIC MFrontGenericInterface
  PUBLIC ${DOLFIN_LIBRARIES})
if(WIN32)
  install(TARGETS MFrontGenericInterface-FEniCS
          DESTINATION bin)
else(WIN32)
  install(TARGETS MFrontGenericInterface-FEniCS
          DESTINATION lib${LIB_SUFFIX})
endif(WIN32)
