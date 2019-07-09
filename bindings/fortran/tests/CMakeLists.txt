STRING(REPLACE "\\;" ";" MGIS_PATH_STRING "$ENV{PATH}")
STRING(REPLACE ";" "\\;" MGIS_PATH_STRING "${MGIS_PATH_STRING}")

add_library(MFrontGenericInterfaceFortranUnitTestingUtilities
  EXCLUDE_FROM_ALL SHARED
  mgis_testing_utilities.f95
  mgis_testing_utilities.c)
target_include_directories(MFrontGenericInterfaceFortranUnitTestingUtilities
 PRIVATE "${CMAKE_BINARY_DIR}/bindings/fortran/modules")
target_link_libraries(MFrontGenericInterfaceFortranUnitTestingUtilities
   PRIVATE MFrontGenericInterface-fortran)

function(test_fortran_bindings test)
  add_executable(${test}_f
    EXCLUDE_FROM_ALL ${test}.f95)
  target_include_directories(${test}_f
    PRIVATE "${CMAKE_BINARY_DIR}/bindings/fortran/modules")
  target_link_libraries(${test}_f
    PRIVATE MFrontGenericInterface-fortran MFrontGenericInterfaceFortranUnitTestingUtilities)
  if(CMAKE_CONFIGURATION_TYPES)
    foreach(conf ${CMAKE_CONFIGURATION_TYPES})
      add_test(NAME ${test}_${conf}_f
	COMMAND ${test}_f CONFIGURATION ${conf})
      set_tests_properties(${test}_${conf}_f
    	PROPERTIES ENVIRONMENT
        "MGIS_TEST_BEHAVIOURS_LIBRARY=$<TARGET_FILE:BehaviourTest>")
    endforeach(conf ${CMAKE_CONFIGURATION_TYPES})
  else(CMAKE_CONFIGURATION_TYPES)
    add_test(NAME ${test}_f COMMAND ${test}_f)
    if((CMAKE_HOST_WIN32) AND (NOT MSYS))
      set_property(TEST ${test}_f
        PROPERTY DEPENDS BehaviourTest
        PROPERTY ENVIRONMENT
	"PATH=$<TARGET_FILE_DIR:MFrontGenericInterface>\;$<TARGET_FILE_DIR:MFrontGenericInterface-c>\;$<TARGET_FILE_DIR:MFrontGenericInterface-fortran>\;${MGIS_PATH_STRING}"
	"MGIS_TEST_BEHAVIOURS_LIBRARY=$<TARGET_FILE:BehaviourTest>")
    else((CMAKE_HOST_WIN32) AND (NOT MSYS))
      set_property(TEST ${test}_f
        PROPERTY DEPENDS BehaviourTest
	PROPERTY ENVIRONMENT
	MGIS_TEST_BEHAVIOURS_LIBRARY="$<TARGET_FILE:BehaviourTest>")
    endif((CMAKE_HOST_WIN32) AND (NOT MSYS))
    add_dependencies(check ${test}_f)
  endif(CMAKE_CONFIGURATION_TYPES)
endfunction(test_fortran_bindings)

test_fortran_bindings(mgis_status)
test_fortran_bindings(MFrontGenericBehaviourInterfaceTest-f)
test_fortran_bindings(MFrontGenericBehaviourInterfaceTest2-f)
test_fortran_bindings(MFrontGenericBehaviourInterfaceTest3-f)
test_fortran_bindings(BehaviourDataTest-f)
test_fortran_bindings(StateTest-f)
test_fortran_bindings(StateTest2-f)
test_fortran_bindings(IntegrateTest-f)
test_fortran_bindings(IntegrateTest2-f)
test_fortran_bindings(IntegrateTest3-f)
test_fortran_bindings(ParameterTest-f)
test_fortran_bindings(BoundsCheckTest-f)
