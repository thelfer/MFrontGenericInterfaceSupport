link_directories("${PROJECT_BINARY_DIR}/bindings/c/src")
link_directories("${PROJECT_BINARY_DIR}/src")

add_executable(MFrontGenericBehaviourInterfaceTest-c
	MFrontGenericBehaviourInterfaceTest-c.cxx)
target_link_libraries(MFrontGenericBehaviourInterfaceTest-c
	PRIVATE MFrontGenericInterface-c)
#add_executable(MFrontGenericBehaviourInterfaceTest2-c
#	MFrontGenericBehaviourInterfaceTest2.cxx)
#target_link_libraries(MFrontGenericBehaviourInterfaceTest2
#	PRIVATE MFrontGenericInterface)

add_test(NAME MFrontGenericBehaviourInterfaceTest-c
 COMMAND MFrontGenericBehaviourInterfaceTest-c
 "$<TARGET_FILE:BehaviourTest>" "gurson")
add_dependencies(check MFrontGenericBehaviourInterfaceTest-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ${test_arg}
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
  set_property(TEST ${test_arg}
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$ENV{PATH}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ${test_arg}
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))
