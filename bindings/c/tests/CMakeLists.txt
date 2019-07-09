STRING(REPLACE "\\;" ";" MGIS_PATH_STRING "$ENV{PATH}")
STRING(REPLACE ";" "\\;" MGIS_PATH_STRING "${MGIS_PATH_STRING}")

add_executable(MFrontGenericBehaviourInterfaceTest-c
  EXCLUDE_FROM_ALL
  MFrontGenericBehaviourInterfaceTest-c.c)
target_compile_definitions(MFrontGenericBehaviourInterfaceTest-c
  PRIVATE -DTFEL_VERSION="${TFEL_VERSION}")
target_link_libraries(MFrontGenericBehaviourInterfaceTest-c
  PRIVATE MFrontGenericInterface-c MFrontGenericInterface)
add_executable(MFrontGenericBehaviourInterfaceTest2-c
  EXCLUDE_FROM_ALL
  MFrontGenericBehaviourInterfaceTest2-c.c)
target_link_libraries(MFrontGenericBehaviourInterfaceTest2-c
  PRIVATE MFrontGenericInterface-c MFrontGenericInterface)
add_executable(MFrontGenericBehaviourInterfaceTest3-c
  EXCLUDE_FROM_ALL
  MFrontGenericBehaviourInterfaceTest3-c.c)
target_link_libraries(MFrontGenericBehaviourInterfaceTest3-c
  PRIVATE MFrontGenericInterface-c MFrontGenericInterface)
add_executable(BehaviourDataTest-c
  EXCLUDE_FROM_ALL
  BehaviourDataTest-c.c)
target_link_libraries(BehaviourDataTest-c
  PRIVATE MFrontGenericInterface-c MFrontGenericInterface)
add_executable(StateTest-c
  EXCLUDE_FROM_ALL StateTest-c.c)
target_link_libraries(StateTest-c
  PRIVATE MFrontGenericInterface-c MFrontGenericInterface)
add_executable(StateTest2-c
  EXCLUDE_FROM_ALL StateTest2-c.c)
target_link_libraries(StateTest2-c
  PRIVATE MFrontGenericInterface-c MFrontGenericInterface)

add_executable(BoundsCheckTest-c
  EXCLUDE_FROM_ALL
  BoundsCheckTest-c.c)
target_compile_definitions(BoundsCheckTest-c
  PRIVATE -DTFEL_VERSION="${TFEL_VERSION}")
target_link_libraries(BoundsCheckTest-c
  PRIVATE MFrontGenericInterface-c MFrontGenericInterface)

add_executable(ParameterTest-c
  EXCLUDE_FROM_ALL
  ParameterTest-c.c)
target_compile_definitions(ParameterTest-c
  PRIVATE -DTFEL_VERSION="${TFEL_VERSION}")
target_link_libraries(ParameterTest-c
  PRIVATE MFrontGenericInterface-c MFrontGenericInterface)

add_executable(IntegrateTest-c
  EXCLUDE_FROM_ALL
  IntegrateTest-c.c)
target_link_libraries(IntegrateTest-c
  PRIVATE MFrontGenericInterface-c MFrontGenericInterface)
add_executable(IntegrateTest2-c
  EXCLUDE_FROM_ALL
  IntegrateTest2-c.c)
target_link_libraries(IntegrateTest2-c
  PRIVATE MFrontGenericInterface-c MFrontGenericInterface)
add_executable(IntegrateTest3-c
  EXCLUDE_FROM_ALL
  IntegrateTest3-c.c)
target_link_libraries(IntegrateTest3-c
  PRIVATE MFrontGenericInterface-c MFrontGenericInterface)

add_test(NAME MFrontGenericBehaviourInterfaceTest-c
 COMMAND MFrontGenericBehaviourInterfaceTest-c
 "$<TARGET_FILE:BehaviourTest>" "Gurson")
add_dependencies(check MFrontGenericBehaviourInterfaceTest-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MFrontGenericBehaviourInterfaceTest-c
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MFrontGenericBehaviourInterfaceTest-c
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME MFrontGenericBehaviourInterfaceTest2-c
 COMMAND MFrontGenericBehaviourInterfaceTest2-c
 "$<TARGET_FILE:BehaviourTest>" "FiniteStrainSingleCrystal")
add_dependencies(check MFrontGenericBehaviourInterfaceTest2-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MFrontGenericBehaviourInterfaceTest2-c
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MFrontGenericBehaviourInterfaceTest2-c
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME MFrontGenericBehaviourInterfaceTest3-c
 COMMAND MFrontGenericBehaviourInterfaceTest3-c
 "$<TARGET_FILE:BehaviourTest>" "FiniteStrainSingleCrystal")
add_dependencies(check MFrontGenericBehaviourInterfaceTest3-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MFrontGenericBehaviourInterfaceTest3-c
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST MFrontGenericBehaviourInterfaceTest3-c
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME BehaviourDataTest-c
 COMMAND BehaviourDataTest-c
 "$<TARGET_FILE:BehaviourTest>" "Gurson")
add_dependencies(check BehaviourDataTest-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST BehaviourDataTest-c
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST BehaviourDataTest-c
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME StateTest-c
 COMMAND StateTest-c
 "$<TARGET_FILE:BehaviourTest>" "Gurson")
add_dependencies(check StateTest-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST StateTest-c
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST StateTest-c
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME StateTest2-c
 COMMAND StateTest2-c
 "$<TARGET_FILE:BehaviourTest>" "Gurson")
add_dependencies(check StateTest2-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST StateTest2-c
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST StateTest2-c
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME BoundsCheckTest-c
 COMMAND BoundsCheckTest-c
 "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check BoundsCheckTest-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST BoundsCheckTest-c
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST BoundsCheckTest-c
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME ParameterTest-c
 COMMAND ParameterTest-c
 "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check ParameterTest-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ParameterTest-c
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST ParameterTest-c
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME IntegrateTest-c
 COMMAND IntegrateTest-c
 "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check IntegrateTest-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest-c
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest-c
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME IntegrateTest2-c
 COMMAND IntegrateTest2-c
 "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check IntegrateTest2-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest2-c
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest2-c
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))

add_test(NAME IntegrateTest3-c
 COMMAND IntegrateTest3-c
 "$<TARGET_FILE:BehaviourTest>")
add_dependencies(check IntegrateTest3-c)
if((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest3-c
    PROPERTY DEPENDS BehaviourTest
    PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;${MGIS_PATH_STRING}")
else((CMAKE_HOST_WIN32) AND (NOT MSYS))
  set_property(TEST IntegrateTest3-c
    PROPERTY DEPENDS BehaviourTest)
endif((CMAKE_HOST_WIN32) AND (NOT MSYS))


