link_directories("${PROJECT_BINARY_DIR}/src"
                 "${TFEL_LIBRARY_PATH}")

mfront_behaviours_check_library(BehaviourTest
  Norton
  Gurson
  FiniteStrainSingleCrystal)

add_executable(MFrontGenericBehaviourInterfaceTest
  EXCLUDE_FROM_ALL
  MFrontGenericBehaviourInterfaceTest.cxx)
target_link_libraries(MFrontGenericBehaviourInterfaceTest
  PRIVATE MFrontGenericInterface)
add_executable(MFrontGenericBehaviourInterfaceTest2
  EXCLUDE_FROM_ALL
  MFrontGenericBehaviourInterfaceTest2.cxx)
target_link_libraries(MFrontGenericBehaviourInterfaceTest2
	PRIVATE MFrontGenericInterface)
add_executable(IntegrateTest
  EXCLUDE_FROM_ALL IntegrateTest.cxx)
target_link_libraries(IntegrateTest
	PRIVATE MFrontGenericInterface)

add_test(NAME MFrontGenericBehaviourInterfaceTest
 COMMAND MFrontGenericBehaviourInterfaceTest
 "$<TARGET_FILE:BehaviourTest>" "Gurson")
add_dependencies(check MFrontGenericBehaviourInterfaceTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MFrontGenericBehaviourInterfaceTest
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$ENV{PATH}")
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
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$ENV{PATH}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MFrontGenericBehaviourInterfaceTest2
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME IntegrateTest
 COMMAND IntegrateTest "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check IntegrateTest)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$ENV{PATH}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))
