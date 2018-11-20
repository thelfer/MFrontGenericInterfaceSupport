subroutine test()
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  ! behaviour
  type(Behaviour) b
  ! hypothesis
  character(len=:), allocatable :: h
  ! source
  character(len=:), allocatable :: s
  ! version
  character(len=:), allocatable :: v
  ! version
  character(len=:), allocatable :: n
  real(kind=8) :: yg_min = 100.d9
  real(kind=8) :: yg_max = 200.d9
  real(kind=8) :: eps = 1.d-14
  real(kind=8) :: upperBound
  real(kind=8) :: lowerPhysicalBound
  integer :: nparams
  logical :: hasBounds, hasLowerBound, hasUpperBound
  logical :: hasPhysicalBounds, hasLowerPhysicalBound, hasUpperPhysicalBound
  logical :: r
  ! start of the check
  call check_status(load_behaviour(b, &
       get_mfront_behaviour_test_library_path(), &
       'BoundsCheckTest', 'Tridimensional'))
  call  check_status(behaviour_get_hypothesis(h, b))
  r = check_string(h, 'Tridimensional', 'invalid hypothesis')
  call check_status(behaviour_get_source(s, b))
  r = check_string(s, 'BoundsCheckTest.mfront', 'invalid source')
  ! version
  call check_status(behaviour_get_tfel_version(v, b))
  r = check_string(v, get_tfel_version(), 'invalid TFEL version')
  ! bounds
  call check_status(behaviour_has_bounds(hasBounds, b, 'YoungModulus'))
  r = check(hasBounds, "'YoungModulus' shall have bounds")
  call check_status(behaviour_has_lower_bound(hasLowerBound, b, "YoungModulus"))
  r = check(hasLowerBound, "'YoungModulus' shall have a lower bound")
  call check_status(behaviour_has_upper_bound(hasUpperBound, b, "YoungModulus"))
  r = check(hasUpperBound, "'YoungModulus' shall have an upper bound")
  call check_status(behaviour_get_upper_bound(upperBound, b, "YoungModulus"))
  write(*,*) 'upperBound: ', upperBound, yg_max, abs(upperBound - yg_max)
  r = check(abs(upperBound - yg_max) < eps * yg_max, &
       "invalid upper bound for 'YoungModulus'")
  ! physical bounds
  call check_status(behaviour_has_physical_bounds(hasPhysicalBounds, b, "YoungModulus"))
  r = check(hasPhysicalBounds, "'YoungModulus' shall have bounds")
  call check_status(behaviour_has_lower_physical_bound(hasLowerPhysicalBound, b, "YoungModulus"))
  r = check(hasLowerPhysicalBound, "'YoungModulus' shall have a lower bound")
  call check_status(behaviour_get_lower_physical_bound(lowerPhysicalBound, &
        b, "YoungModulus"))
  r = check(abs(lowerPhysicalBound) < eps * yg_min, &
       "invalid physical lower bound for 'YoungModulus'")
  call check_status(behaviour_has_upper_physical_bound(hasUpperPhysicalBound, b, "YoungModulus"))
  r = check(.not. hasUpperPhysicalBound, "'YoungModulus' shall have an upper bound")
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

