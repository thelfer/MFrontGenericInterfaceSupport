subroutine test()
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  type(Behaviour) b
  type(BehaviourData) d
  type(State) s0,s1
  real(kind=8) de ! strain increment
  real(kind=8) p0 ! equivalent viscoplastic strain value at the
                  ! beginning of the time step
  real(kind=8) p1 ! equivalent viscoplastic strain value at the
                  ! end of the time step
  real(kind=8), dimension(21) :: p_ref ! referecence values of the
                                       ! equivalent plastic strain
  real(kind=8), dimension(21) ::p ! computed values of the equivalent
                                  ! plastic strain
  real(kind=8), dimension(6) :: e
  integer i ! index
  integer o ! offset of the equivalent viscoplastic strain
  integer ri ! returned value of a behaviour integration
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
  call check_status(allocate_behaviour_data(d,b))
  ! state at the beginning of the time step
  call check_status(behaviour_data_get_state_0(s0, d))
  ! state at the end of the time step
  call check_status(behaviour_data_get_state_1(s1, d))
  ! setting the temperature at the beginning of the time step
  call check_status(state_set_external_state_variable_by_name( &
       s0, 'Temperature', 293.15d0))
  ! getting the offset of the equivalent plastic strain
  call check_status(behaviour_get_internal_state_variable_offset( &
       o, b, 'EquivalentViscoplasticStrain'))
  ! setting the time increment
  call check_status(behaviour_data_set_time_increment(d, 180d0))
  ! getting the addresses where the equivalent plastic strain is stored
  call check_status(state_get_internal_state_variable_by_offset(p0, s0, o))
  call check_status(state_get_internal_state_variable_by_offset(p1, s1, o))
  ! strain at the end of the time step
  e = (/ de, 0d0, 0d0, 0d0, 0d0, 0d0 /)
  call check_status(state_set_gradient_by_name(s1, 'Strain', e))
  p(1) = p0
  do i = 1, 20
     call check_status(integrate(ri, d, b))
     r = check(ri.eq.1, 'integration failed')
     call check_status(update_behaviour_data(d))
     e(1) = e(1) + de
     call check_status(state_set_gradient_by_name(s1, 'Strain', e))
     call check_status(state_get_internal_state_variable_by_offset(p1, s1, o))
     p(i + 1) = p1
  end do
  do i = 1, 21
     r = check(abs(p(i)-p_ref(i)) < 1.d-12, &
          'invalid value for the equivalent plastic strain ' // &
          ' vs ' )
  end do
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

