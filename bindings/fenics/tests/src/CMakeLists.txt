add_custom_command(OUTPUT MGISSmallStrainFormulation3D.cpp MGISSmallStrainFormulation3D.h
                   COMMAND ffc -l dolfin -O -fno-evaluate_basis -fsplit ${CMAKE_CURRENT_SOURCE_DIR}/MGISSmallStrainFormulation3D.ufl
                   MAIN_DEPENDENCY MGISSmallStrainFormulation3D.ufl)
add_custom_command(OUTPUT MGISSmallStrainFormulation2D.cpp MGISSmallStrainFormulation2D.h
                   COMMAND ffc -l dolfin -O -fno-evaluate_basis -fsplit ${CMAKE_CURRENT_SOURCE_DIR}/MGISSmallStrainFormulation2D.ufl
                   MAIN_DEPENDENCY MGISSmallStrainFormulation2D.ufl)

add_library(MFrontGenericInterfaceFEniCSTestingUtilities SHARED
  EXCLUDE_FROM_ALL
  FEniCSTestingUtilities.cxx
  MGISSmallStrainFormulation2D.cpp
  MGISSmallStrainFormulation3D.cpp)
if (COMPILER_DEFAULT_VISIBILITY_FLAG)
  target_compile_options(MFrontGenericInterfaceFEniCSTestingUtilities
    PRIVATE "${COMPILER_DEFAULT_VISIBILITY_FLAG}")
endif(COMPILER_DEFAULT_VISIBILITY_FLAG)
target_include_directories(MFrontGenericInterfaceFEniCSTestingUtilities
  PUBLIC "${CMAKE_SOURCE_DIR}/include" "${CMAKE_SOURCE_DIR}/bindings/fenics/tests/include")
target_include_directories(MFrontGenericInterfaceFEniCSTestingUtilities
  SYSTEM PUBLIC ${DOLFIN_INCLUDE_DIRS} ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})
target_link_libraries(MFrontGenericInterfaceFEniCSTestingUtilities
  PUBLIC ${DOLFIN_LIBRARIES})

function(test_fenics_bindings test)
  add_executable(${test}_fenics
    EXCLUDE_FROM_ALL ${test}.cxx)
  target_include_directories(${test}_fenics
    SYSTEM PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
  target_include_directories(${test}_fenics
    PRIVATE "${CMAKE_SOURCE_DIR}/bindings/fenics/tests/include")
  target_link_libraries(${test}_fenics
    PRIVATE MFrontGenericInterface-FEniCS
    PRIVATE MFrontGenericInterfaceFEniCSTestingUtilities)
  if(CMAKE_CONFIGURATION_TYPES)
    foreach(conf ${CMAKE_CONFIGURATION_TYPES})
      add_test(NAME ${test}_${conf}_f
	COMMAND ${test}_fenics CONFIGURATION ${conf})
      set_tests_properties(${test}_${conf}_f
    	PROPERTIES ENVIRONMENT
        MGIS_FENICS_TEST_DATA="${CMAKE_CURRENT_BINARY_DIR}/data"
	MGIS_TEST_BEHAVIOURS_LIBRARY="$<TARGET_FILE:BehaviourTest>"
	MGIS_TEST_TFEL_VERSION=${TFEL_VERSION})
    endforeach(conf ${CMAKE_CONFIGURATION_TYPES})
  else(CMAKE_CONFIGURATION_TYPES)
    add_test(NAME ${test}_fenics COMMAND ${test}_fenics)
    set_property(TEST ${test}_fenics
      PROPERTY ENVIRONMENT
      MGIS_FENICS_TEST_DATA=${CMAKE_SOURCE_DIR}/bindings/fenics/tests/data
      MGIS_TEST_BEHAVIOURS_LIBRARY=$<TARGET_FILE:BehaviourTest>
      MGIS_TEST_TFEL_VERSION=${TFEL_VERSION})
    add_dependencies(check ${test}_fenics)
  endif(CMAKE_CONFIGURATION_TYPES)
endfunction(test_fenics_bindings)

# test_fenics_bindings(ElasticityUniaxialTensileTestImposedStrainPlaneStress-exx)

test_fenics_bindings(ElasticityUniaxialTensileTestImposedStress3D-sxx)
test_fenics_bindings(ElasticityUniaxialTensileTestImposedStrain3D-exx)
test_fenics_bindings(ElasticityUniaxialTensileTestImposedStrain3D-eyy)
test_fenics_bindings(ElasticityUniaxialTensileTestImposedStrain3D-ezz)
test_fenics_bindings(ElasticityUniaxialTensileTestImposedStrain3D-exz)
test_fenics_bindings(ElasticityUniaxialTensileTestImposedStrain3D-exy)
test_fenics_bindings(ElasticityUniaxialTensileTestImposedStrain3D-eyz)
test_fenics_bindings(PlasticityUniaxialTensileTestImposedStress3D-sxx)
test_fenics_bindings(StandardElastoPlasticityPlasticityTest11-cyclic_E)
# test_fenics_bindings(PlasticCylinderExpansion)
