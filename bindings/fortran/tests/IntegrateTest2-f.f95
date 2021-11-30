subroutine test()
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  type(Behaviour) b
  type(MaterialDataManager) m
  type(MaterialStateManager) s0,s1
  real(kind=8) de ! strain increment
  real(kind=8) p0 ! equivalent viscoplastic strain value at the
                  ! beginning of the time step
  real(kind=8) p1 ! equivalent viscoplastic strain value at the
                  ! end of the time step
  real(kind=8) dt ! time step
  real(kind=8), dimension(21) :: p_ref ! referecence values of the
                                       ! equivalent plastic strain
  real(kind=8), dimension(21) ::pi ! computed values of the equivalent
                                   ! plastic strain
  real(kind=8), dimension(21) ::pe ! computed values of the equivalent
                                   ! plastic strain
  real(kind=8), dimension(6) :: e
  real(kind=8), dimension (:,:), pointer :: g => null()
  real(kind=8), dimension (:,:), pointer :: isvs => null()
  integer :: n_g     ! gradients stride
  integer :: n_isvs  ! gradients stride
  integer :: n = 100 ! index
  integer :: i, j    ! index
  integer :: o       ! offset of the equivalent viscoplastic strain
  integer :: ri      ! returned value of a behaviour integration
  logical :: r
  ! reference values of the equivalent plastic strain
  p_ref = (/ 0d0, 1.3523277308229d-11, &
       1.0955374667213d-07, &
       5.5890770166084d-06, &
       3.2392193670428d-05, &
       6.645865307584d-05, &
       9.9676622883138d-05, &
       0.00013302758358953d0, &
       0.00016635821069889d0, &
       0.00019969195920296d0, &
       0.00023302522883648d0, &
       0.00026635857194317d0, &
       0.000299691903777d0, &
       0.0003330252373404d0, &
       0.00036635857063843d0, &
       0.00039969190397718d0, &
       0.00043302523730968d0, &
       0.00046635857064314d0, &
       0.00049969190397646d0, &
       0.00053302523730979d0, &
       0.00056635857064313d0 /)
  ! strain increment
  de = 5.d-5
  ! start of the check
  call check_status(load_behaviour(b, &
       get_mfront_behaviour_test_library_path(), &
       'Norton', 'Tridimensional'))
  ! allocating data
  call check_status(create_material_data_manager(m, b, n))
  ! state at the beginning of the time step
  call check_status(material_data_manager_get_state_0(s0, m))
  ! state at the end of the time step
  call check_status(material_data_manager_get_state_1(s1, m))
  ! 
  call check_status(behaviour_get_internal_state_variable_offset_by_name( &
       o, b, 'EquivalentViscoplasticStrain'))
  ! initialize the external state variable
  call check_status(material_state_manager_set_us_external_state_variable( &
       s1, "Temperature", 293.15d0))
  ! copy s1 in s0
  call check_status(update_material_data_manager(m))
  ! Getting a pointer to the gradients.
  call check_status(material_state_manager_get_gradients(g,s1))
  call check_status(material_state_manager_get_gradients_stride(n_g,s1))
  r = check(n_g == 6, 'invalid gradient stride')
  !  Getting a pointer to the internal state variables.
  call check_status(material_state_manager_get_internal_state_variables(isvs,s1))
  call check_status(material_state_manager_get_internal_state_variables_stride(n_isvs,s1))
  r = check(n_isvs == 7, 'invalid gradient stride')
  ! initializing the strain
  do i = 1, N
     g(:,i) = (/ de, 0d0, 0d0, 0d0, 0d0, 0d0 /)
  end do  
  ! initial value of the equivalent plastic strain
  pi(1) = isvs(o, 1)
  pe(1) = isvs(o, n)
  ! time step
  dt = 180d0
  ! integration
  do j = 1, 20
     call check_status(integrate_material_data_manager_part(ri, m, INTEGRATION_NO_TANGENT_OPERATOR, dt, 1, n))
     r = check(ri.eq.1, 'integration failed')
     call check_status(update_material_data_manager(m))
     ! updating the strain
     do i = 1, N
        g(1,i) = g(1,i) + de
     end do
     pi(j+1) = isvs(o, 1)
     pe(j+1) = isvs(o, n)
  end do
  do i = 1, 21
     r = check(abs(pi(i)-p_ref(i)) < 1.d-12, &
          'invalid value for the equivalent plastic strain ')
     r = check(abs(pi(i)-p_ref(i)) < 1.d-12, &
          'invalid value for the equivalent plastic strain ')
  end do
  ! free ressources
  call check_status(free_material_data_manager(m))
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

