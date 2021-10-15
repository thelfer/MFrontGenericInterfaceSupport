subroutine test()
  use mgis
  use mgis_behaviour
  use mgis_model
  use mgis_testing_utilities
  implicit none
  type(ThreadPool) p
  type(Behaviour) mo
  type(MaterialDataManager) m
  type(MaterialStateManager) s0,s1
  real(kind=8) A
  real(kind=8) eps
  real(kind=8) t ! current time
  real(kind=8) dt ! time step
  real(kind=8) x_ref
  real(kind=8), dimension(11) ::xi ! computed values of x
  real(kind=8), dimension(11) ::xe ! computed values of x
  real(kind=8), dimension (:,:), pointer :: isvs => null()
  integer :: n_g     ! gradients stride
  integer :: n_isvs  ! gradients stride
  integer :: n = 100 ! index
  integer :: i, j    ! index
  integer :: o       ! offset of the x
  integer :: ri      ! returned value of a behaviour integration
  logical :: r
  !
  eps = 1e-10
  ! creation of the thread pool
  call check_status(create_thread_pool(p, 2))
  ! start of the check
  call check_status(load_behaviour(mo, &
       get_mfront_model_test_library_path(), &
       'ode_rk54', 'Tridimensional'))
  call check_status(behaviour_get_parameter_default_value(A, mo, 'A'))
  ! allocating data
  call check_status(create_material_data_manager(m, mo, n))
  ! state at the beginning of the time step
  call check_status(material_data_manager_get_state_0(s0, m))
  ! state at the end of the time step
  call check_status(material_data_manager_get_state_1(s1, m))
  ! 
  call check_status(behaviour_get_internal_state_variable_offset( &
       o, mo, 'x'))
  ! initialize the external state variable
  call check_status(material_state_manager_set_uniform_external_state_variable( &
       s1, "Temperature", 293.15d0))
  !  Getting a pointer to the internal state variables.
  call check_status(material_state_manager_get_internal_state_variables(isvs,s1))
  ! initial value of the x
  do i = 1, N
     isvs(o, 1) = 1
  end do  
  ! copy s1 in s0
  call check_status(update_material_data_manager(m))
  ! initial value of the x
  xi(1) = isvs(o, 1)
  xe(1) = isvs(o, n)
  ! time step
  dt = 0.1d0
  ! integration
  do j = 1, 10
     call check_status(integrate_material_data_manager( &
          ri, p, m, INTEGRATION_NO_TANGENT_OPERATOR, dt))
     r = check(ri.eq.1, 'integration failed')
     call check_status(update_material_data_manager(m))
     xi(j+1) = isvs(o, 1)
     xe(j+1) = isvs(o, n)
  end do
  ! check results
  t = 0d0
  do i = 1, 11
     x_ref = exp(-A * t)
     r = check(abs(xi(i)-x_ref) < eps, &
          'invalid value for x')
     r = check(abs(xi(i)-x_ref) < eps, &
          'invalid value for x')
     t = t + dt
  end do
  ! free ressources
  call check_status(free_material_data_manager(m))
  call check_status(free_behaviour(mo))
  call check_status(free_thread_pool(p))
end subroutine test

program main
  use mgis_testing_utilities
  call test()
  call tests_summary()
  if(.not. status) then
     stop -1
  end if
end program main

