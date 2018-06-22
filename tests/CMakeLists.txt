link_directories("${PROJECT_BINARY_DIR}/src")
add_library(BehaviourTest SHARED
  Gurson.cxx
  FiniteStrainSingleCrystal.cxx)
add_executable(MFrontGenericBehaviourInterfaceTest
	MFrontGenericBehaviourInterfaceTest.cxx)
target_link_libraries(MFrontGenericBehaviourInterfaceTest
	PRIVATE MFrontGenericInterface)
add_executable(MFrontGenericBehaviourInterfaceTest2
	MFrontGenericBehaviourInterfaceTest2.cxx)
target_link_libraries(MFrontGenericBehaviourInterfaceTest2
	PRIVATE MFrontGenericInterface)

add_test(NAME MFrontGenericBehaviourInterfaceTest
 COMMAND MFrontGenericBehaviourInterfaceTest
 "$<TARGET_FILE:BehaviourTest>" "gurson")
add_dependencies(check MFrontGenericBehaviourInterfaceTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MFrontGenericBehaviourInterfaceTest
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$ENV{PATH}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ${test_arg}
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME MFrontGenericBehaviourInterfaceTest2
 COMMAND MFrontGenericBehaviourInterfaceTest2
 "$<TARGET_FILE:BehaviourTest>" "FiniteStrainSingleCrystal")
add_dependencies(check MFrontGenericBehaviourInterfaceTest2)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MFrontGenericBehaviourInterfaceTest2
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$ENV{PATH}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ${test_arg}
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))


