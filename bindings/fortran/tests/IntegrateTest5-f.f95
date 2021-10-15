subroutine test()
  use mgis
  use mgis_behaviour
  use mgis_model
  use mgis_testing_utilities
  implicit none
  type(Behaviour) mo
  type(BehaviourData) d
  type(State) s1
  real(kind=8) A
  real(kind=8) eps
  real(kind=8) t ! current time
  real(kind=8) dt ! time step
  real(kind=8) x_ref
  real(kind=8), dimension(11) ::xvalues ! computed values of x
  real(kind=8), dimension (:), pointer :: isvs => null()
  integer :: i, j    ! index
  integer :: o       ! offset of the x
  integer :: ri      ! returned value of a behaviour integration
  logical :: r
  !
  eps = 1e-10
  ! start of the check
  call check_status(load_behaviour(mo, &
       get_mfront_model_test_library_path(), &
       'ode_rk54', 'Tridimensional'))
  call check_status(behaviour_get_parameter_default_value(A, mo, 'A'))
  ! allocating data
  call check_status(allocate_behaviour_data(d, mo))
  ! state at the end of the time step
  call check_status(behaviour_data_get_state_1(s1, d))
  ! 
  call check_status(behaviour_get_internal_state_variable_offset( &
       o, mo, 'x'))
  ! initialize the external state variable
  call check_status(state_set_external_state_variable_by_name( &
       s1, "Temperature", 293.15d0))
  !  Getting the value of the internal state variable.
  call check_status(state_get_internal_state_variables(isvs, s1, mo))
  ! initial value of the x
  isvs(o) = 1
  ! copy s1 in s0
  call check_status(update_behaviour_data(d))
  ! initial value of the x
  xvalues(1) = isvs(o)
  ! time step
  dt = 0.1d0
  ! setting the time increment
  call check_status(behaviour_data_set_time_increment(d, dt))
  ! integration
  do j = 1, 10
     call check_status(integrate(ri, d, mo))
     r = check(ri.eq.1, 'integration failed')
     call check_status(update_behaviour_data(d))
     !  Getting the value of the internal state variable.
     xvalues(j+1) = isvs(o)
  end do
  ! check results
  t = 0d0
  do i = 1, 11
     x_ref = exp(-A * t)
     r = check(abs(xvalues(i)-x_ref) < eps, &
          'invalid value for x')
     t = t + dt
  end do
  ! free ressources
  call check_status(free_behaviour_data(d))
  call check_status(free_behaviour(mo))
end subroutine test

program main
  use mgis_testing_utilities
  call test()
  call tests_summary()
  if(.not. status) then
     stop -1
  end if
end program main

