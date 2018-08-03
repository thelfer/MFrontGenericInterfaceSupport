link_directories("${PROJECT_BINARY_DIR}/bindings/c/src")
link_directories("${PROJECT_BINARY_DIR}/src")

add_executable(MFrontGenericBehaviourInterfaceTest-c
	MFrontGenericBehaviourInterfaceTest-c.c)
target_link_libraries(MFrontGenericBehaviourInterfaceTest-c
	PRIVATE MFrontGenericInterface-c)
add_executable(MFrontGenericBehaviourInterfaceTest2-c
	MFrontGenericBehaviourInterfaceTest2-c.c)
target_link_libraries(MFrontGenericBehaviourInterfaceTest2-c
	PRIVATE MFrontGenericInterface-c)
add_executable(BehaviourDataTest-c
	BehaviourDataTest-c.c)
target_link_libraries(BehaviourDataTest-c
	PRIVATE MFrontGenericInterface-c)
add_executable(StateTest-c
	StateTest-c.c)
target_link_libraries(StateTest-c
	PRIVATE MFrontGenericInterface-c)

add_test(NAME MFrontGenericBehaviourInterfaceTest-c
 COMMAND MFrontGenericBehaviourInterfaceTest-c
 "$<TARGET_FILE:BehaviourTest>" "gurson")
add_dependencies(check MFrontGenericBehaviourInterfaceTest-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MFrontGenericBehaviourInterfaceTest-c
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;$ENV{PATH}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ${test_arg}
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME MFrontGenericBehaviourInterfaceTest2-c
 COMMAND MFrontGenericBehaviourInterfaceTest2-c
 "$<TARGET_FILE:BehaviourTest>" "FiniteStrainSingleCrystal")
add_dependencies(check MFrontGenericBehaviourInterfaceTest2-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MFrontGenericBehaviourInterfaceTest2-c
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;$ENV{PATH}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ${test_arg}
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME BehaviourDataTest-c
 COMMAND BehaviourDataTest-c
 "$<TARGET_FILE:BehaviourTest>" "gurson")
add_dependencies(check BehaviourDataTest-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST BehaviourDataTest-c
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;$ENV{PATH}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ${test_arg}
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME StateTest-c
 COMMAND StateTest-c
 "$<TARGET_FILE:BehaviourTest>" "gurson")
add_dependencies(check StateTest-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST StateTest-c
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;$ENV{PATH}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ${test_arg}
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

