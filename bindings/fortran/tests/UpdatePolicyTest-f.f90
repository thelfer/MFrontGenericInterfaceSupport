! revert copies the state at the beginning of the time step in the state
! at the end of the time step, except the values set with NOUPDATE.
! The Poisson ratio is null, so the behaviour computes sig_xx = E eto_xx.
subroutine test(p, Eref)
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  integer, intent(in) :: p
  real(kind=8), intent(in) :: Eref
  type(Behaviour) b
  type(MaterialDataManager) m
  type(MaterialStateManager) s0, s1
  real(kind=8), dimension(:,:), pointer :: g => null()
  real(kind=8), dimension(:,:), pointer :: sig => null()
  real(kind=8) :: e = 1d-3
  integer ri
  logical :: r
  call check_status(load_behaviour(b, &
       get_mfront_behaviour_test_library_path(), &
       'PostProcessingTest', 'Tridimensional'))
  call check_status(create_material_data_manager(m, b, 1))
  call check_status(material_data_manager_get_state_0(s0, m))
  call check_status(material_data_manager_get_state_1(s1, m))
  call check_status(material_state_manager_set_uniform_scalar_material_property( &
       s0, 'YoungModulus', 100d9))
  call check_status(material_state_manager_set_uniform_scalar_material_property( &
       s0, 'PoissonRatio', 0d0))
  call check_status(material_state_manager_set_us_external_state_variable( &
       s0, 'Temperature', 293.15d0))
  call check_status(material_state_manager_set_uniform_scalar_material_property2( &
       s1, 'YoungModulus', 200d9, p))
  call check_status(revert_material_data_manager(m))
  call check_status(material_state_manager_get_gradients(g, s1))
  g(:, 1) = (/ e, 0d0, 0d0, 0d0, 0d0, 0d0 /)
  call check_status(integrate_material_data_manager_part( &
       ri, m, INTEGRATION_NO_TANGENT_OPERATOR, 1d0, 1, 1))
  r = check(ri == 1, 'integration failed')
  call check_status(material_state_manager_get_thermodynamic_forces(sig, s1))
  r = check(abs(sig(1, 1) - Eref * e) < 1d-8 * Eref * e, &
       'invalid stress (update policy not respected)')
  ! free ressources
  call check_status(free_material_data_manager(m))
  call check_status(free_behaviour(b))
end subroutine test

program main
  use mgis_behaviour
  use mgis_testing_utilities
  call test(NOUPDATE, 200d9)
  call test(UPDATE, 100d9)
  call tests_summary()
  if(.not. status) then
     stop -1
  end if
end program main
