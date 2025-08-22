STRING(REPLACE "\\;" ";" MGIS_PATH_STRING "$ENV{PATH}")
STRING(REPLACE ";" "\\;" MGIS_PATH_STRING "${MGIS_PATH_STRING}")

mfront_material_properties_check_library(MaterialPropertyTest
  MaterialPropertyBoundCheck
  MaterialPropertyBoundCheck2
  MaterialPropertyBoundCheck3
  MaterialPropertyErrnoCheck
  MaterialPropertyErrnoCheck2
  Inconel600_YoungModulus)

mfront_behaviours_check_library(BehaviourTest
  Elasticity
  OrthotropicElasticity
  Norton
  NortonFailure
  Plasticity
  Gurson
  BoundsCheckTest
  ParameterTest
  FiniteStrainSingleCrystal
  StandardElastoViscoPlasticityPlasticityTest11
  TensorialExternalStateVariableTest
  InitializeFunctionTest
  PostProcessingTest
  ComputeSpeedOfSound)

mfront_behaviours_check_library(ModelTest
  ode_rk54)

#

add_executable(IntegrationFailureDebugTest
  EXCLUDE_FROM_ALL
  IntegrationFailureDebugTest.cxx)
target_link_libraries(IntegrationFailureDebugTest
  PRIVATE MFrontGenericInterface)

add_executable(LoadMaterialPropertyTest
  EXCLUDE_FROM_ALL
  LoadMaterialPropertyTest.cxx)
target_link_libraries(LoadMaterialPropertyTest
  PRIVATE MFrontGenericInterface)
target_compile_definitions(LoadMaterialPropertyTest
  PRIVATE -DTFEL_VERSION="${TFEL_VERSION}")

add_test(NAME LoadMaterialPropertyTest
 COMMAND LoadMaterialPropertyTest
 "$<TARGET_FILE:MaterialPropertyTest>")
add_dependencies(check LoadMaterialPropertyTest)

# Test on behaviours

add_executable(MFrontGenericBehaviourInterfaceTest
  EXCLUDE_FROM_ALL
  MFrontGenericBehaviourInterfaceTest.cxx)
target_link_libraries(MFrontGenericBehaviourInterfaceTest
  PRIVATE MFrontGenericInterface)
target_compile_definitions(MFrontGenericBehaviourInterfaceTest
  PRIVATE -DTFEL_VERSION="${TFEL_VERSION}")

add_executable(MFrontGenericBehaviourInterfaceTest2
  EXCLUDE_FROM_ALL
  MFrontGenericBehaviourInterfaceTest2.cxx)
target_compile_definitions(MFrontGenericBehaviourInterfaceTest2
  PRIVATE -DTFEL_VERSION="${TFEL_VERSION}")
target_link_libraries(MFrontGenericBehaviourInterfaceTest2
  PRIVATE MFrontGenericInterface)

add_executable(MFrontGenericBehaviourInterfaceTest3
  EXCLUDE_FROM_ALL
  MFrontGenericBehaviourInterfaceTest3.cxx)
target_compile_definitions(MFrontGenericBehaviourInterfaceTest3
  PRIVATE -DTFEL_VERSION="${TFEL_VERSION}")
target_link_libraries(MFrontGenericBehaviourInterfaceTest3
  PRIVATE MFrontGenericInterface)

add_executable(BoundsCheckTest
  EXCLUDE_FROM_ALL
  BoundsCheckTest.cxx)
target_link_libraries(BoundsCheckTest
  PRIVATE MFrontGenericInterface)
target_compile_definitions(BoundsCheckTest
  PRIVATE -DTFEL_VERSION="${TFEL_VERSION}")

add_executable(ParameterTest
  EXCLUDE_FROM_ALL
  ParameterTest.cxx)
target_link_libraries(ParameterTest
  PRIVATE MFrontGenericInterface)
target_compile_definitions(ParameterTest
  PRIVATE -DTFEL_VERSION="${TFEL_VERSION}")

add_executable(IntegrateTest
  EXCLUDE_FROM_ALL IntegrateTest.cxx)
target_link_libraries(IntegrateTest
	PRIVATE MFrontGenericInterface)
add_executable(IntegrateTest2
  EXCLUDE_FROM_ALL IntegrateTest2.cxx)
target_link_libraries(IntegrateTest2
	PRIVATE MFrontGenericInterface)
add_executable(IntegrateTest2b
  EXCLUDE_FROM_ALL IntegrateTest2b.cxx)
target_link_libraries(IntegrateTest2b
	PRIVATE MFrontGenericInterface)
add_executable(IntegrateTest3
  EXCLUDE_FROM_ALL IntegrateTest3.cxx)
target_link_libraries(IntegrateTest3
	PRIVATE MFrontGenericInterface)
add_executable(IntegrateTest3b
  EXCLUDE_FROM_ALL IntegrateTest3b.cxx)
target_link_libraries(IntegrateTest3b
	PRIVATE MFrontGenericInterface)
add_executable(IntegrateTest4
  EXCLUDE_FROM_ALL IntegrateTest4.cxx)
target_link_libraries(IntegrateTest4
	PRIVATE MFrontGenericInterface)
add_executable(IntegrateTest5
  EXCLUDE_FROM_ALL IntegrateTest5.cxx)
target_link_libraries(IntegrateTest5
	PRIVATE MFrontGenericInterface)
add_executable(IntegrateDebugTest
  EXCLUDE_FROM_ALL IntegrateDebugTest.cxx)
target_link_libraries(IntegrateDebugTest
	PRIVATE MFrontGenericInterface)
add_executable(RotateFunctionsTest
  EXCLUDE_FROM_ALL RotateFunctionsTest.cxx)
target_link_libraries(RotateFunctionsTest
	PRIVATE MFrontGenericInterface)

add_executable(ExternalStateVariableTest
  EXCLUDE_FROM_ALL ExternalStateVariableTest.cxx)
target_link_libraries(ExternalStateVariableTest
	PRIVATE MFrontGenericInterface)
target_compile_definitions(ExternalStateVariableTest
  PRIVATE -DTFEL_VERSION="${TFEL_VERSION}")

add_executable(InitializeFunctionTest
  EXCLUDE_FROM_ALL InitializeFunctionTest.cxx)
target_link_libraries(InitializeFunctionTest
	PRIVATE MFrontGenericInterface)
target_compile_definitions(InitializeFunctionTest
  PRIVATE -DTFEL_VERSION="${TFEL_VERSION}")

add_executable(PostProcessingTest
  EXCLUDE_FROM_ALL PostProcessingTest.cxx)
target_link_libraries(PostProcessingTest
	PRIVATE MFrontGenericInterface)
target_compile_definitions(PostProcessingTest
  PRIVATE -DTFEL_VERSION="${TFEL_VERSION}")

add_executable(ComputeSpeedOfSoundTest
  EXCLUDE_FROM_ALL ComputeSpeedOfSoundTest.cxx)
target_link_libraries(ComputeSpeedOfSoundTest
	PRIVATE MFrontGenericInterface)

add_executable(ComputeSpeedOfSoundTest2
  EXCLUDE_FROM_ALL ComputeSpeedOfSoundTest2.cxx)
target_link_libraries(ComputeSpeedOfSoundTest2
	PRIVATE MFrontGenericInterface)

add_test(NAME MFrontGenericBehaviourInterfaceTest
 COMMAND MFrontGenericBehaviourInterfaceTest
 "$<TARGET_FILE:BehaviourTest>" "Gurson")
add_dependencies(check MFrontGenericBehaviourInterfaceTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MFrontGenericBehaviourInterfaceTest
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MFrontGenericBehaviourInterfaceTest
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME MFrontGenericBehaviourInterfaceTest2
 COMMAND MFrontGenericBehaviourInterfaceTest2
 "$<TARGET_FILE:BehaviourTest>" "FiniteStrainSingleCrystal")
add_dependencies(check MFrontGenericBehaviourInterfaceTest2)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MFrontGenericBehaviourInterfaceTest2
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MFrontGenericBehaviourInterfaceTest2
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME MFrontGenericBehaviourInterfaceTest3
 COMMAND MFrontGenericBehaviourInterfaceTest3
 "$<TARGET_FILE:BehaviourTest>" "FiniteStrainSingleCrystal")
add_dependencies(check MFrontGenericBehaviourInterfaceTest3)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MFrontGenericBehaviourInterfaceTest3
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MFrontGenericBehaviourInterfaceTest3
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME BoundsCheckTest
 COMMAND BoundsCheckTest
 "$<TARGET_FILE:BehaviourTest>" "BoundsCheckTest")
add_dependencies(check BoundsCheckTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST BoundsCheckTest
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST BoundsCheckTest
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME ParameterTest
 COMMAND ParameterTest
 "$<TARGET_FILE:BehaviourTest>" "ParameterTest")
add_dependencies(check ParameterTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ParameterTest
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ParameterTest
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME IntegrateTest
 COMMAND IntegrateTest "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check IntegrateTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME IntegrateTest2
 COMMAND IntegrateTest2 "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check IntegrateTest2)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest2
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest2
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME IntegrateTest2b
 COMMAND IntegrateTest2b "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check IntegrateTest2b)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest2b
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest2b
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME IntegrateTest3
 COMMAND IntegrateTest3 "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check IntegrateTest3)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest3
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest3
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME IntegrateTest3b
 COMMAND IntegrateTest3b "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check IntegrateTest3b)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest3b
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest3b
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME RotateFunctionsTest
 COMMAND RotateFunctionsTest "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check RotateFunctionsTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST RotateFunctionsTest
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST RotateFunctionsTest
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME IntegrateTest4
 COMMAND IntegrateTest4 "$<TARGET_FILE:ModelTest>")
add_dependencies(check IntegrateTest4)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest4
    PROPERTY DEPENDS ModelTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest4
    PROPERTY DEPENDS ModelTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME IntegrateTest5
 COMMAND IntegrateTest5 "$<TARGET_FILE:ModelTest>")
add_dependencies(check IntegrateTest5)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest5
    PROPERTY DEPENDS ModelTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest5
    PROPERTY DEPENDS ModelTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME IntegrationFailureDebugTest
 COMMAND IntegrationFailureDebugTest "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check IntegrationFailureDebugTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrationFailureDebugTest
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrationFailureDebugTest
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME IntegrateDebugTest
 COMMAND IntegrateDebugTest "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check IntegrateDebugTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateDebugTest
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateDebugTest
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME ExternalStateVariableTest
  COMMAND ExternalStateVariableTest "$<TARGET_FILE:BehaviourTest>"  "TensorialExternalStateVariableTest")
add_dependencies(check ExternalStateVariableTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ExternalStateVariableTest
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ExternalStateVariableTest
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME InitializeFunctionTest
  COMMAND InitializeFunctionTest "$<TARGET_FILE:BehaviourTest>"  "InitializeFunctionTest")
add_dependencies(check InitializeFunctionTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST InitializeFunctionTest
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST InitializeFunctionTest
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME PostProcessingTest
  COMMAND PostProcessingTest "$<TARGET_FILE:BehaviourTest>"  "PostProcessingTest")
add_dependencies(check PostProcessingTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST PostProcessingTest
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST PostProcessingTest
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME ComputeSpeedOfSoundTest
 COMMAND ComputeSpeedOfSoundTest "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check ComputeSpeedOfSoundTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ComputeSpeedOfSoundTest
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ComputeSpeedOfSoundTest
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME ComputeSpeedOfSoundTest2
 COMMAND ComputeSpeedOfSoundTest2 "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check ComputeSpeedOfSoundTest2)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ComputeSpeedOfSoundTest2
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ComputeSpeedOfSoundTest2
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))
