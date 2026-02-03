STRING(REPLACE "\\;" ";" MGIS_PATH_STRING "$ENV{PATH}")
STRING(REPLACE ";" "\\;" MGIS_PATH_STRING "${MGIS_PATH_STRING}")

if(MGIS_HAVE_MFRONT_SUPPORT)

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
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST LoadMaterialPropertyTest
    PROPERTY DEPENDS MaterialPropertyTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST LoadMaterialPropertyTest
    PROPERTY DEPENDS MaterialPropertyTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

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

add_executable(UpdatePolicyTest
  EXCLUDE_FROM_ALL UpdatePolicyTest.cxx)
target_link_libraries(UpdatePolicyTest
	PRIVATE MFrontGenericInterface)

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

add_test(NAME UpdatePolicyTest
 COMMAND UpdatePolicyTest "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check UpdatePolicyTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST UpdatePolicyTest
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST UpdatePolicyTest
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

endif(MGIS_HAVE_MFRONT_SUPPORT)

add_executable(InvalidResultTest
 EXCLUDE_FROM_ALL
 InvalidResultTest.cxx)
target_link_libraries(InvalidResultTest
  PRIVATE MFrontGenericInterface)
add_test(NAME InvalidResultTest COMMAND InvalidResultTest)
add_dependencies(check InvalidResultTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST InvalidResultTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_executable(ErrorBacktraceTest
 EXCLUDE_FROM_ALL
 ErrorBacktraceTest.cxx)
target_link_libraries(ErrorBacktraceTest
  PRIVATE MFrontGenericInterface)
add_test(NAME ErrorBacktraceTest COMMAND ErrorBacktraceTest)
add_dependencies(check ErrorBacktraceTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ErrorBacktraceTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

if(enable-mgis-function)

if(MGIS_HAVE_TFEL)

add_executable(SharedSpaceTest
 EXCLUDE_FROM_ALL
 SharedSpaceTest.cxx)
target_link_libraries(SharedSpaceTest
  PRIVATE MFrontGenericInterface tfel::TFELTests)
add_test(NAME SharedSpaceTest COMMAND SharedSpaceTest)
add_dependencies(check SharedSpaceTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST SharedSpaceTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_executable(FunctionTest
 EXCLUDE_FROM_ALL
 FunctionTest.cxx)
target_link_libraries(FunctionTest
  PRIVATE MFrontGenericInterface tfel::TFELTests)
add_test(NAME FunctionTest COMMAND FunctionTest)
add_dependencies(check FunctionTest)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.2.0)
    target_compile_definitions(FunctionTest
       PRIVATE MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS)
  endif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.2.0)
endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
if(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    target_compile_definitions(FunctionTest
       PRIVATE MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS)
endif()
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST FunctionTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_executable(EvaluatorsTest
 EXCLUDE_FROM_ALL
 EvaluatorsTest.cxx)
target_link_libraries(EvaluatorsTest
  PRIVATE MFrontGenericInterface tfel::TFELTests)
add_test(NAME EvaluatorsTest COMMAND EvaluatorsTest)
add_dependencies(check EvaluatorsTest)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.2.0)
    target_compile_definitions(EvaluatorsTest
       PRIVATE MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS)
  endif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.2.0)
endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
if(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    target_compile_definitions(EvaluatorsTest
       PRIVATE MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS)
endif()
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST EvaluatorsTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_executable(UniformEvaluatorsTest
 EXCLUDE_FROM_ALL
 UniformEvaluatorsTest.cxx)
target_link_libraries(UniformEvaluatorsTest
  PRIVATE MFrontGenericInterface tfel::TFELTests)
add_test(NAME UniformEvaluatorsTest COMMAND UniformEvaluatorsTest)
add_dependencies(check UniformEvaluatorsTest)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.2.0)
    target_compile_definitions(UniformEvaluatorsTest
       PRIVATE MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS)
  endif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.2.0)
endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
if(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    target_compile_definitions(UniformEvaluatorsTest
       PRIVATE MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS)
endif()
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST UniformEvaluatorsTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_executable(CoalescedMemoryAccessFunctionViewTest
 EXCLUDE_FROM_ALL
 CoalescedMemoryAccessFunctionViewTest.cxx)
target_link_libraries(CoalescedMemoryAccessFunctionViewTest
  PRIVATE MFrontGenericInterface tfel::TFELTests)
add_test(NAME CoalescedMemoryAccessFunctionViewTest
         COMMAND CoalescedMemoryAccessFunctionViewTest)
add_dependencies(check CoalescedMemoryAccessFunctionViewTest)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.2.0)
    target_compile_definitions(CoalescedMemoryAccessFunctionViewTest
       PRIVATE MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS)
  endif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.2.0)
endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
if(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    target_compile_definitions(CoalescedMemoryAccessFunctionViewTest
       PRIVATE MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS)
endif()
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST CoalescedMemoryAccessFunctionViewTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_executable(StridedCoalescedMemoryAccessFunctionViewTest
 EXCLUDE_FROM_ALL
 StridedCoalescedMemoryAccessFunctionViewTest.cxx)
target_link_libraries(StridedCoalescedMemoryAccessFunctionViewTest
  PRIVATE MFrontGenericInterface tfel::TFELTests)
add_test(NAME StridedCoalescedMemoryAccessFunctionViewTest
         COMMAND StridedCoalescedMemoryAccessFunctionViewTest)
add_dependencies(check StridedCoalescedMemoryAccessFunctionViewTest)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.2.0)
    target_compile_definitions(StridedCoalescedMemoryAccessFunctionViewTest
       PRIVATE MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS)
  endif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.2.0)
endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
if(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    target_compile_definitions(StridedCoalescedMemoryAccessFunctionViewTest
       PRIVATE MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS)
endif()
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST StridedCoalescedMemoryAccessFunctionViewTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_executable(MechanicalEvaluatorsTest
 EXCLUDE_FROM_ALL
 MechanicalEvaluatorsTest.cxx)
target_link_libraries(MechanicalEvaluatorsTest
  PRIVATE MFrontGenericInterface tfel::TFELTests)
add_test(NAME MechanicalEvaluatorsTest COMMAND MechanicalEvaluatorsTest)
add_dependencies(check MechanicalEvaluatorsTest)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.2.0)
    target_compile_definitions(MechanicalEvaluatorsTest
       PRIVATE MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS)
  endif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.2.0)
endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
if(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    target_compile_definitions(MechanicalEvaluatorsTest
       PRIVATE MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS)
endif()
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MechanicalEvaluatorsTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))


add_executable(TensorialFunctionTest
 EXCLUDE_FROM_ALL
 TensorialFunctionTest.cxx)
target_link_libraries(TensorialFunctionTest
  PRIVATE MFrontGenericInterface tfel::TFELTests)
add_test(NAME TensorialFunctionTest COMMAND TensorialFunctionTest)
add_dependencies(check TensorialFunctionTest)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.2.0)
    target_compile_definitions(TensorialFunctionTest
       PRIVATE MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS)
  endif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.2.0)
endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
if(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    target_compile_definitions(TensorialFunctionTest
       PRIVATE MGIS_DISABLE_CONSTEXPR_FUNCTION_TESTS)
endif()
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST TensorialFunctionTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_executable(ReduceTest
 EXCLUDE_FROM_ALL
 ReduceTest.cxx)
target_link_libraries(ReduceTest
  PRIVATE MFrontGenericInterface tfel::TFELTests)
add_test(NAME ReduceTest COMMAND ReduceTest)
add_dependencies(check ReduceTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ReduceTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_executable(MaterialFunctionManagerTest
  EXCLUDE_FROM_ALL MaterialFunctionManagerTest.cxx)
target_link_libraries(MaterialFunctionManagerTest
	PRIVATE MFrontGenericInterface tfel::TFELTests)

if(MGIS_HAVE_MFRONT_SUPPORT)

add_test(NAME MaterialFunctionManagerTest
 COMMAND MaterialFunctionManagerTest "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check MaterialFunctionManagerTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MaterialFunctionManagerTest
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MaterialFunctionManagerTest
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

endif(MGIS_HAVE_MFRONT_SUPPORT)

endif(MGIS_HAVE_TFEL)

if(MGIS_HDF5_SUPPORT)

  add_executable(HDF5Test EXCLUDE_FROM_ALL HDF5Test.cxx)
  target_link_libraries(HDF5Test PRIVATE MFrontGenericInterface)
  add_test(NAME HDF5Test
   COMMAND HDF5Test "$<TARGET_FILE:BehaviourTest>")
  add_dependencies(check HDF5Test)
  if((CMAKE_HOST_WIN32) AND (NOT MSYS))
    set_property(TEST HDF5Test
      PROPERTY DEPENDS BehaviourTest
      PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;${MGIS_PATH_STRING}")
  else((CMAKE_HOST_WIN32) AND (NOT MSYS))
    set_property(TEST HDF5Test
      PROPERTY DEPENDS BehaviourTest)
  endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

endif(MGIS_HDF5_SUPPORT)

endif(enable-mgis-function)
