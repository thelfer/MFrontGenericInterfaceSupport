subroutine test()
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  type(Behaviour) b
  type(BehaviourData) d
  type(State) s0,s1
  real(kind=8) T ! temperature value
  logical :: r
  ! start of the check
  call check_status(load_behaviour(b, &
       get_mfront_behaviour_test_library_path(), &
       'Gurson', 'Tridimensional'))
  call check_status(allocate_behaviour_data(d,b))
  ! state at the beginning of the time step
  call check_status(behaviour_data_get_state_0(s0, d));
  ! state at the end of the time step
  call check_status(behaviour_data_get_state_1(s1, d));
  ! setting the temperature at the beginning of the time step
  call check_status(state_set_external_state_variable_by_name( &
       s0, "Temperature", 293.15d0));
  ! s0 is copied in s1
  call check_status(revert_behaviour_data(d));
  call check_status(state_get_external_state_variable_by_name(T, s0, &
       "Temperature"));
  r =check(abs(T - 293.15d0) < 1e-8,"invalid temperature value");
  call check_status(state_get_external_state_variable_by_name(T, s1, &
       "Temperature"));
  r =check(abs(T - 293.15d0) < 1e-8,"invalid temperature value");
  ! free ressources
  call check_status(free_behaviour_data(d))
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

