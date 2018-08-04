link_directories("${PROJECT_BINARY_DIR}/bindings/c/src")
link_directories("${PROJECT_BINARY_DIR}/src")

add_executable(MFrontGenericBehaviourInterfaceTest-c
  EXCLUDE_FROM_ALL
  MFrontGenericBehaviourInterfaceTest-c.c)
target_link_libraries(MFrontGenericBehaviourInterfaceTest-c
  PRIVATE MFrontGenericInterface-c)
add_executable(MFrontGenericBehaviourInterfaceTest2-c
  EXCLUDE_FROM_ALL
  MFrontGenericBehaviourInterfaceTest2-c.c)
target_link_libraries(MFrontGenericBehaviourInterfaceTest2-c
  PRIVATE MFrontGenericInterface-c)
add_executable(BehaviourDataTest-c
  EXCLUDE_FROM_ALL
  BehaviourDataTest-c.c)
target_link_libraries(BehaviourDataTest-c
  PRIVATE MFrontGenericInterface-c)
add_executable(StateTest-c
  EXCLUDE_FROM_ALL StateTest-c.c)
target_link_libraries(StateTest-c
  PRIVATE MFrontGenericInterface-c)
add_executable(StateTest2-c
  EXCLUDE_FROM_ALL StateTest2-c.c)
target_link_libraries(StateTest2-c
  PRIVATE MFrontGenericInterface-c)
add_executable(IntegrateTest-c
  EXCLUDE_FROM_ALL
  IntegrateTest-c.c)
target_link_libraries(IntegrateTest-c
  PRIVATE MFrontGenericInterface-c)

add_test(NAME MFrontGenericBehaviourInterfaceTest-c
 COMMAND MFrontGenericBehaviourInterfaceTest-c
 "$<TARGET_FILE:BehaviourTest>" "Gurson")
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
 "$<TARGET_FILE:BehaviourTest>" "Gurson")
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
 "$<TARGET_FILE:BehaviourTest>" "Gurson")
add_dependencies(check StateTest-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST StateTest-c
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;$ENV{PATH}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ${test_arg}
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME StateTest2-c
 COMMAND StateTest2-c
 "$<TARGET_FILE:BehaviourTest>" "Gurson")
add_dependencies(check StateTest2-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST StateTest2-c
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;$ENV{PATH}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ${test_arg}
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME IntegrateTest-c
 COMMAND IntegrateTest-c
 "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check IntegrateTest-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest-c
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;$ENV{PATH}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ${test_arg}
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

