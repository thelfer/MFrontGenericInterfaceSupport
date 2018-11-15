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
  ! hypothesis
  character(len=:), allocatable :: h
  ! source
  character(len=:), allocatable :: s
  ! version
  character(len=:), allocatable :: v
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
  ! hypothesis
  call  check_status(behaviour_get_hypothesis(h, b))
  r = check_string(h, 'Tridimensional', 'invalid hypothesis')
  ! source
  call check_status(behaviour_get_source(s, b))
  r = check_string(s, 'Gurson.mfront', 'invalid source')
  ! version
  call check_status(behaviour_get_tfel_version(v, b));
  r = check_string(v, get_tfel_version(), "invalid TFEL version");
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
