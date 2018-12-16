add_custom_command(OUTPUT Plas3D.cpp Plas3D.h
                   COMMAND ffc -l dolfin -O -fno-evaluate_basis -fsplit ${CMAKE_CURRENT_SOURCE_DIR}/Plas3D.ufl
                   MAIN_DEPENDENCY Plas3D.ufl)

function(test_fenics_bindings test)
  add_executable(${test}_fenics
    EXCLUDE_FROM_ALL ${test}.cxx Plas3D.cpp)
  target_include_directories(${test}_fenics
   SYSTEM PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
  target_link_libraries(${test}_fenics
    PRIVATE MFrontGenericInterface-FEniCS)
  if(CMAKE_CONFIGURATION_TYPES)
    foreach(conf ${CMAKE_CONFIGURATION_TYPES})
      add_test(NAME ${test}_${conf}_f
	COMMAND ${test}_fenics CONFIGURATION ${conf})
      set_tests_properties(${test}_${conf}_f
    	PROPERTIES ENVIRONMENT
        MGIS_TEST_BEHAVIOURS_LIBRARY="$<TARGET_FILE:BehaviourTest>"
	MGIS_TEST_TFEL_VERSION=${TFEL_VERSION})
    endforeach(conf ${CMAKE_CONFIGURATION_TYPES})
  else(CMAKE_CONFIGURATION_TYPES)
    add_test(NAME ${test}_fenics COMMAND ${test}_fenics)
    set_property(TEST ${test}_fenics
      PROPERTY ENVIRONMENT
      MGIS_TEST_BEHAVIOURS_LIBRARY=$<TARGET_FILE:BehaviourTest>
      MGIS_TEST_TFEL_VERSION=${TFEL_VERSION})
    add_dependencies(check ${test}_fenics)
  endif(CMAKE_CONFIGURATION_TYPES)
endfunction(test_fenics_bindings)

test_fenics_bindings(uniaxial-tensile-test)
