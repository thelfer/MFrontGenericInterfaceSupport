subroutine test()
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  type(behaviour) b
  logical :: r
  ! library
  character(len=:), allocatable :: l
  ! behaviour name
  character(len=:), allocatable :: bn
  ! behaviour name
  character(len=:), allocatable :: fn
  ! hypothesis
  character(len=:), allocatable :: h
  ! source
  character(len=:), allocatable :: s
  ! version
  character(len=:), allocatable :: v
  ! name of the first internal state variable
  character(len=:), allocatable :: eel
  ! name of the second internal state variable
  character(len=:), allocatable :: p
  ! name of the third internal state variable
  character(len=:), allocatable :: pm
  ! name of the fourth internal state variable
  character(len=:), allocatable :: f
  ! name of the first external state variable
  character(len=:), allocatable :: T
  ! number of material properties
  integer mps_size
  ! number of material properties
  integer isvs_size
  ! number of external state variables
  integer esvs_size
  ! type of the first internal state variable
  integer eel_t
  ! type of the second internal state variable
  integer p_t
  ! type of the third internal state variable
  integer pm_t
  ! type of the fourth internal state variable
  integer f_t
  ! type of the first external state variable
  integer T_t
  ! start of the check
  call check_status(load_behaviour(b, &
       get_mfront_behaviour_test_library_path(), &
       'Gurson', 'Tridimensional'))
  call check_status(behaviour_get_library(l, b))
  r = check_string(l, get_mfront_behaviour_test_library_path(), &
       'invalid library')
  ! behaviour name
  call check_status(behaviour_get_behaviour_name(bn, b))
  r = check_string(bn, 'Gurson', 'invalid behaviour name')
  ! function name
  call check_status(behaviour_get_function_name(fn, b))
  r = check_string(fn, 'Gurson_Tridimensional', 'invalid function name')
  ! hypothesis
  call  check_status(behaviour_get_hypothesis(h, b))
  r = check_string(h, 'Tridimensional', 'invalid hypothesis')
  ! source
  call check_status(behaviour_get_source(s, b))
  r = check_string(s, 'Gurson.mfront', 'invalid source')
  ! version
  call check_status(behaviour_get_tfel_version(v, b));
  r = check_string(v, get_tfel_version(), "invalid TFEL version")
  ! material properties
  call check_status(behaviour_get_number_of_material_properties(mps_size, b))
  r = check(mps_size .eq. 0, "invalid number of material properties")
  ! internal state variable
  call check_status(behaviour_get_number_of_internal_state_variables(isvs_size, b))
  if (check(isvs_size .eq. 4, "invalid number of internal state variable")) then
     call check_status(behaviour_get_internal_state_variable_name(eel, b, 1))
     r = check_string(eel, "ElasticStrain", "invalid internal state variable name")
     call check_status(behaviour_get_internal_state_variable_type(eel_t, b, 1))
     r = check(eel_t == STENSOR,"invalid type for internal state variable 'ElasticStrain'")
     call check_status(behaviour_get_internal_state_variable_name(p, b, 2))
     r = check_string(p, "EquivalentPlasticStrain", "invalid internal state variable name")
     call check_status(behaviour_get_internal_state_variable_type(p_t, b, 2))
     r = check(p_t == SCALAR,"invalid type for internal state variable 'EquivalentPlasticStrain'")
     call check_status(behaviour_get_internal_state_variable_name(pm, b, 3))
     r = check_string(pm, "MatrixEquivalentPlasticStrain", "invalid internal state variable name")
     call check_status(behaviour_get_internal_state_variable_type(pm_t, b, 3))
     r = check(pm_t == SCALAR,"invalid type for internal state variable 'MatrixEquivalentPlasticStrain'")
     call check_status(behaviour_get_internal_state_variable_name(f, b, 4))
     r = check_string(f, "Porosity", "invalid internal state variable name")
     call check_status(behaviour_get_internal_state_variable_type(f_t, b, 4))
     r = check(f_t == SCALAR,"invalid type for internal state variable 'Porosity'")
  end if
  ! external state variable
  call check_status(behaviour_get_number_of_external_state_variables(esvs_size, b))
  if (check(esvs_size .eq. 1, "invalid number of external state variable")) then
     call check_status(behaviour_get_external_state_variable_name(T, b, 1))
     r = check_string(T, "Temperature", "invalid external state variable name")
     call check_status(behaviour_get_external_state_variable_type(T_t, b, 1))
     r = check(T_t == SCALAR,"invalid type for external state variable 'Temperature'")
  end if
  ! free ressources
  call check_status(free_behaviour(b))
end subroutine test

program main
  use mgis_testing_utilities
  call test()
  call tests_summary()
  if(.not. status) then
     stop -1
  end if
end program main
