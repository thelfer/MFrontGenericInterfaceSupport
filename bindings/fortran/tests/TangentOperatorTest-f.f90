! tests of the tangent and prediction operators returned by the integrations
! of a material data manager and of a behaviour data, with the `TangentKind`
! behaviour whose operator identifies the type of integration. The Fortran
! arrays of the tangent operators are the transposes of the operators.
subroutine test()
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  type(ThreadPool) p
  type(Behaviour) b
  type(MaterialDataManager) m
  type(MaterialStateManager) s0, s1
  type(BehaviourData) d
  type(State) ds0, ds1
  real(kind=8), dimension(:,:), pointer :: g => null()
  real(kind=8), dimension(:,:), pointer :: isvs => null()
  real(kind=8), dimension(:,:,:), pointer :: K => null()
  real(kind=8), dimension(:,:), pointer :: Kd => null()
  real(kind=8), dimension(:), pointer :: Ka => null()
  real(kind=8), parameter :: kv = 10 ! value of the material property k
  integer, dimension(7) :: types
  integer :: it, mode, idx, i, j, t
  integer :: ri ! returned value of a behaviour integration
  logical :: r, treated
  types = (/ PREDICTION_TANGENT_OPERATOR, PREDICTION_SECANT_OPERATOR, &
       PREDICTION_ELASTIC_OPERATOR, INTEGRATION_ELASTIC_OPERATOR, &
       INTEGRATION_SECANT_OPERATOR, INTEGRATION_TANGENT_OPERATOR, &
       INTEGRATION_CONSISTENT_TANGENT_OPERATOR /)
  call check_status(create_thread_pool(p, 2))
  call check_status(load_behaviour(b, &
       get_mfront_behaviour_test_library_path(), &
       'TangentKind', 'Tridimensional'))
  do it = 1, 7
     t = types(it)
     ! mode 1: points 2 and 3, mode 2: thread pool on all the points
     do mode = 1, 2
        call check_status(create_material_data_manager(m, b, 4))
        call check_status(material_data_manager_get_state_0(s0, m))
        call check_status(material_data_manager_get_state_1(s1, m))
        call check_status(material_state_manager_set_uniform_scalar_material_property( &
             s0, 'k', kv))
        call check_status(material_state_manager_set_uniform_scalar_material_property( &
             s1, 'k', kv))
        call check_status(material_state_manager_set_us_external_state_variable( &
             s0, 'Temperature', 293.15d0))
        call check_status(material_state_manager_set_us_external_state_variable( &
             s1, 'Temperature', 293.15d0))
        call check_status(material_state_manager_get_gradients(g, s1))
        do idx = 1, 4
           g(1, idx) = 1.d-3 * idx
        end do
        if (mode == 1) then
           call check_status(integrate_material_data_manager_part( &
                ri, m, t, 1.d0, 2, 3))
        else
           call check_status(integrate_material_data_manager(ri, p, m, t, 1.d0))
        end if
        r = check(ri == 1, 'integration failed')
        call check_status(material_data_manager_get_tangent_operator(K, m))
        call check_status(material_state_manager_get_internal_state_variables(isvs, s1))
        do idx = 1, 4
           ! the points outside the range are left untouched
           treated = (mode == 2) .or. ((idx >= 2) .and. (idx <= 3))
           do i = 1, 6
              do j = 1, 6
                 if (treated) then
                    r = check(abs(K(j, i, idx) - expected_operator(t, kv, i, j)) &
                         < 1.d-12 * kv, 'invalid value of the operator')
                 else
                    r = check(abs(K(j, i, idx)) < 1.d-12 * kv, &
                         'invalid value of the operator')
                 end if
              end do
           end do
           ! a prediction leaves the internal state variables untouched
           if (treated .and. (t > 0)) then
              r = check(abs(isvs(1, idx) - 1) < 1.d-14, &
                   'invalid value of the internal state variable')
           else
              r = check(abs(isvs(1, idx)) < 1.d-14, &
                   'invalid value of the internal state variable')
           end if
        end do
        call check_status(free_material_data_manager(m))
     end do
     ! behaviour data
     call check_status(allocate_behaviour_data(d, b))
     call check_status(behaviour_data_get_state_0(ds0, d))
     call check_status(behaviour_data_get_state_1(ds1, d))
     call check_status(state_set_material_property_by_name(ds0, 'k', 7.d0))
     call check_status(state_set_material_property_by_name(ds1, 'k', 7.d0))
     call check_status(state_set_scalar_external_state_variable_by_name( &
          ds0, 'Temperature', 293.15d0))
     call check_status(state_set_scalar_external_state_variable_by_name( &
          ds1, 'Temperature', 293.15d0))
     call check_status(behaviour_data_set_time_increment(d, 1.d0))
     call check_status(behaviour_data_get_tangent_operator_array(Ka, d))
     Ka(1) = t
     call check_status(integrate(ri, d, b))
     r = check(ri == 1, 'integration failed')
     call check_status(behaviour_data_get_tangent_operator(Kd, d))
     do i = 1, 6
        do j = 1, 6
           r = check(abs(Kd(j, i) - expected_operator(t, 7.d0, i, j)) < 1.d-12 * 7, &
                'invalid value of the operator of the behaviour data')
        end do
     end do
     call check_status(free_behaviour_data(d))
  end do
  call check_status(free_behaviour(b))
  call check_status(free_thread_pool(p))
contains
  ! component (i, j) of the operator computed by the `TangentKind` behaviour
  ! for the type of integration t and the value k of the material property
  function expected_operator(t, k, i, j) result(v)
    implicit none
    integer, intent(in) :: t, i, j
    real(kind=8), intent(in) :: k
    real(kind=8) :: v
    v = 0
    if (i == j) then
       v = t * k
    end if
    if ((t > 0) .and. (i == 1) .and. (j == 2)) then
       v = v + k / 2
    end if
  end function expected_operator
end subroutine test

program main
  use mgis_testing_utilities
  call test()
  call tests_summary()
  if(.not. status) then
     stop -1
  end if
end program main
