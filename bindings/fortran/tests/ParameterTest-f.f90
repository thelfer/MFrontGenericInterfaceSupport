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
  real(kind=8) :: pv
  real(kind=8) :: yg = 150d9
  real(kind=8) :: nu = 0.3d0
  real(kind=8) :: eps = 1.d-14
  integer :: nparams
  logical :: r
  ! start of the check
  call check_status(load_behaviour(b, &
       get_mfront_behaviour_test_library_path(), &
       'ParameterTest', 'Tridimensional'))
  call  check_status(behaviour_get_hypothesis(h, b))
  r = check_string(h, 'Tridimensional', 'invalid hypothesis')
  call check_status(behaviour_get_source(s, b))
  r = check_string(s, 'ParameterTest.mfront', 'invalid source')
  ! version
  call check_status(behaviour_get_tfel_version(v, b))
  r = check_string(v, get_tfel_version(), "invalid TFEL version")
  call check_status(behaviour_get_number_of_parameters(nparams, b))
  if (check(nparams == 6, 'invalid number of parameters')) then
     call check_status(behaviour_get_parameter_name(n, b, 1))
     r = check_string(n, 'YoungModulus', 'invalid first parameter')
     call check_status(behaviour_get_parameter_default_value(pv, b, n))
     r = check(abs(pv-yg)<eps*yg,'invalid "YoungModulus" default value')
     call check_status(behaviour_get_parameter_name(n, b, 2))
     r = check_string(n, 'PoissonRatio', 'invalid second parameter')
     call check_status(behaviour_get_parameter_default_value(pv, b, n))
     r = check(abs(pv-nu)<eps*nu,"invalid 'PoissonRatio' default value")
     call check_status(behaviour_get_parameter_name(n, b, 3))
     r = check_string(n, 'ParametersArray[0]', 'invalid third parameter')
     call check_status(behaviour_get_parameter_name(n, b, 4))
     r = check_string(n, 'ParametersArray[1]', 'invalid fourth parameter')
     call check_status(behaviour_get_parameter_name(n, b, 5))
     r = check_string(n, 'minimal_time_step_scaling_factor', 'invalid fifth parameter')
     call check_status(behaviour_get_parameter_name(n, b, 6))
     r = check_string(n, 'maximal_time_step_scaling_factor', 'invalid sixth parameter')
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

