add_test(NAME IntegrateTest_jl
  COMMAND ${Julia_EXECUTABLE} --startup-file=no
          "${CMAKE_CURRENT_SOURCE_DIR}/IntegrateTest.jl")
set_property(TEST IntegrateTest_jl
  PROPERTY ENVIRONMENT
  MGIS_JULIA_LIBRARY=$<TARGET_FILE:mgis-julia>
  MGIS_TEST_BEHAVIOURS_LIBRARY=$<TARGET_FILE:BehaviourTest>)
