subroutine check_external_state_variable(b, expected_name, expected_type, expected_offset, position)
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  type(Behaviour) b
  character(*) :: expected_name
  integer :: expected_type
  integer :: expected_offset
  integer :: position
  character(len=:), allocatable :: name
  integer :: type
  integer :: offset
  logical :: r
  call check_status(behaviour_get_external_state_variable_name(name, b, position))
  r =  check_string(name, expected_name, "invalid external state variable name")
  call check_status(behaviour_get_external_state_variable_type(type, b, position))
  r = check(type == expected_type, "invalid external state variable type")
  call check_status(behaviour_get_external_state_variable_offset(offset, b, position))
  write(*,*) offset, expected_offset
  r = check(offset == expected_offset, "invalid external state variable offset")
end subroutine check_external_state_variable

subroutine check_behaviour(b)
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  type(Behaviour) b
  character(len=:), allocatable :: h
  character(len=:), allocatable :: s
  character(len=:), allocatable :: v
  logical :: r
  integer :: n
  !
  call  check_status(behaviour_get_hypothesis(h, b))
  r = check_string(h, 'Tridimensional', 'invalid hypothesis')
  call check_status(behaviour_get_source(s, b))
  r = check_string(s, 'TensorialExternalStateVariableTest.mfront', &
       'invalid source')
  ! version
  call check_status(behaviour_get_tfel_version(v, b))
  r = check_string(v, get_tfel_version(), "invalid TFEL version")
  ! external state variables
  call check_status(behaviour_get_number_of_external_state_variables(n, b));
  r = check(n == 10, "invalid number of external state variables")
  ! type of the external state variables
  call check_external_state_variable(b, "Temperature", SCALAR, 1, 1)
  call check_external_state_variable(b, "v_esv", VECTOR, 2, 2)
  call check_external_state_variable(b, "v2_esv[0]", VECTOR, 5, 3)
  call check_external_state_variable(b, "v2_esv[1]", VECTOR, 8, 4)
end subroutine check_behaviour

subroutine check_behaviour_data(b)
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  ! behaviour
  type(Behaviour) b
  type(BehaviourData) d
  type(State) s1
  real(kind=8), dimension(3) ::v_esv_values
  logical :: r
  !
  v_esv_values = (/ 1, 2, 3 /)
  ! start of the check
  call check_status(allocate_behaviour_data(d, b))
  !
  call check_status(behaviour_data_get_state_1(s1, d))
  !  call check_status(state_set_external_state_variable_by_name(s1, "v_esv", v_esv_values));

  ! free ressources
  call check_status(free_behaviour_data(d))
endsubroutine check_behaviour_data

program main
  use mgis_behaviour
  use mgis_testing_utilities
  ! behaviour
  type(Behaviour) b
  ! start of the check
  call check_status(load_behaviour(b, &
       get_mfront_behaviour_test_library_path(), &
       'TensorialExternalStateVariableTest', 'Tridimensional'))
  call check_behaviour(b)
  call tests_summary()
  ! free ressources
  call check_status(free_behaviour(b))
  if(.not. status) then
     stop -1
  end if
end program main

