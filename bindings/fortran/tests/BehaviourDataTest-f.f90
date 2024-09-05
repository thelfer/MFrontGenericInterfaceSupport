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

