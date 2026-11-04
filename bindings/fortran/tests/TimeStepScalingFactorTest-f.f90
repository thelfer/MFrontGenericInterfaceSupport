subroutine test()
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  type(ThreadPool) p
  type(Behaviour) b
  type(MaterialDataManager) m
  type(MaterialStateManager) s1
  type(BehaviourIntegrationOptions) o
  real(kind=8) :: tsf ! time step scaling factor proposed by the behaviour
  integer :: ri       ! returned value of a behaviour integration
  logical :: r
  call check_status(create_thread_pool(p, 2))
  call check_status(load_behaviour(b, &
       get_mfront_behaviour_test_library_path(), &
       'TimeStepScalingFactorTest', 'Tridimensional'))
  call check_status(create_material_data_manager(m, b, 3))
  call check_status(material_data_manager_get_state_1(s1, m))
  call check_status(material_state_manager_set_uniform_scalar_material_property( &
       s1, "ProposedTimeStepScalingFactor", 0.5d0))
  call check_status(material_state_manager_set_us_external_state_variable( &
       s1, "Temperature", 293.15d0))
  call check_status(update_material_data_manager(m))
  call check_status(create_behaviour_integration_options(o))
  call check_status(behaviour_integration_options_set_integration_type( &
       o, INTEGRATION_NO_TANGENT_OPERATOR))
  ! all the integration points
  call check_status(integrate_material_data_manager_with_options_2( &
       ri, tsf, p, m, o, 1.d0))
  r = check(abs(tsf - 0.5d0) < 1.d-14, 'invalid time step scaling factor')
  call check_status(integrate_debug_material_data_manager_with_options_2( &
       ri, tsf, p, m, o, 1.d0))
  r = check(abs(tsf - 0.5d0) < 1.d-14, 'invalid time step scaling factor')
  ! last integration point only
  call check_status(integrate_material_data_manager_part_2( &
       ri, tsf, m, INTEGRATION_NO_TANGENT_OPERATOR, 1.d0, 3, 3))
  r = check(abs(tsf - 0.5d0) < 1.d-14, 'invalid time step scaling factor')
  call check_status(integrate_debug_material_data_manager_part_2( &
       ri, tsf, m, INTEGRATION_NO_TANGENT_OPERATOR, 1.d0, 3, 3))
  r = check(abs(tsf - 0.5d0) < 1.d-14, 'invalid time step scaling factor')
  call check_status(free_behaviour_integration_options(o))
  call check_status(free_material_data_manager(m))
  call check_status(free_behaviour(b))
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
