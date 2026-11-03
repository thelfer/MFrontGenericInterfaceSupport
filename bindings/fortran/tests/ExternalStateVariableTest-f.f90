subroutine check_external_state_variable(b, expected_name, expected_type, expected_offset, position)
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  type(Behaviour) b
  character(*) :: expected_name
  integer :: expected_type
  integer :: expected_offset
  integer :: position
  character(len=:), allocatable :: name
  integer :: type
  integer :: offset
  logical :: r
  call check_status(behaviour_get_external_state_variable_name(name, b, position))
  r =  check_string(name, expected_name, "invalid external state variable name")
  call check_status(behaviour_get_external_state_variable_type(type, b, position))
  r = check(type == expected_type, "invalid external state variable type")
  call check_status(behaviour_get_external_state_variable_offset(offset, b, position))
  r = check(offset == expected_offset, "invalid external state variable offset")
end subroutine check_external_state_variable

subroutine check_behaviour(b)
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  type(Behaviour) b
  character(len=:), allocatable :: h
  character(len=:), allocatable :: s
  character(len=:), allocatable :: v
  logical :: r
  integer :: n
  !
  call  check_status(behaviour_get_hypothesis(h, b))
  r = check_string(h, 'Tridimensional', 'invalid hypothesis')
  call check_status(behaviour_get_source(s, b))
  r = check_string(s, 'TensorialExternalStateVariableTest.mfront', &
       'invalid source')
  ! version
  call check_status(behaviour_get_tfel_version(v, b))
  r = check_string(v, get_tfel_version(), "invalid TFEL version")
  ! external state variables
  call check_status(behaviour_get_number_of_external_state_variables(n, b));
  r = check(n == 10, "invalid number of external state variables")
  ! type of the external state variables
  call check_external_state_variable(b, "Temperature", SCALAR, 1, 1)
  call check_external_state_variable(b, "v_esv", VECTOR, 2, 2)
  call check_external_state_variable(b, "v2_esv[0]", VECTOR, 5, 3)
  call check_external_state_variable(b, "v2_esv[1]", VECTOR, 8, 4)
end subroutine check_behaviour

subroutine check_behaviour_data(b)
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  ! behaviour
  type(Behaviour) b
  type(BehaviourData) d
  type(State) s1
  real(kind=8), dimension(3) ::v_esv_values
  real(kind=8), dimension(:), pointer :: isvs
  real(kind=8), dimension(:), pointer :: esv
  real(kind=8) :: x
  type(mgis_status) :: st
  integer :: i
  integer :: ri
  logical :: r
  !
  v_esv_values = (/ 1, 2, 3 /)
  ! start of the check
  call check_status(allocate_behaviour_data(d, b))
  !
  call check_status(behaviour_data_get_state_1(s1, d))
  call check_status(state_set_external_state_variable_by_name(s1, "v_esv", v_esv_values))
  call check_status(state_get_external_state_variable_values_by_name(esv, s1, "v_esv"))
  r = check(size(esv) == 3, "invalid size of the external state variable")
  do i = 1, 3
     r = check(abs(esv(i) - v_esv_values(i)) < 1.d-14, &
          "invalid external state variable value")
  end do
  st = state_get_external_state_variable_by_name(x, s1, "v_esv")
  r = check(st % exit_status .ne. MGIS_SUCCESS, &
       "the scalar getter shall refuse a vector")
  call check_status(integrate(ri, d, b))
  ! the behaviour copies the external state variables in the internal ones
  call check_status(state_get_internal_state_variables(isvs, s1, b))
  do i = 1, 3
     r = check(abs(isvs(i) - v_esv_values(i)) < 1.d-14, &
          "invalid internal state variable value")
  end do
  call check_status(state_set_internal_state_variable_by_name( &
       s1, "v_isv", (/ 4d0, 5d0, 6d0 /)))
  do i = 1, 3
     r = check(abs(isvs(i) - (3 + i)) < 1.d-14, &
          "invalid internal state variable value")
  end do
  ! free ressources
  call check_status(free_behaviour_data(d))
endsubroutine check_behaviour_data

subroutine check_material_data_manager(b)
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  type(Behaviour) b
  type(MaterialDataManager) m
  type(MaterialStateManager) s0, s1
  ! a Young modulus and values of v_esv per integration point
  real(kind=8), dimension(2) :: E = (/ 100d9, 200d9 /)
  real(kind=8), dimension(6) :: v_esv = (/ 1d0, 2d0, 3d0, 4d0, 5d0, 6d0 /)
  real(kind=8), dimension(9) :: zeros = 0d0
  real(kind=8), dimension(:,:), pointer :: g => null()
  real(kind=8), dimension(:,:), pointer :: sig => null()
  real(kind=8), dimension(:,:), pointer :: isvs => null()
  character(len=:), allocatable :: h
  character(len=:), allocatable :: n
  type(mgis_status) :: st
  integer :: i, j, t, vs, nesv, ri
  logical :: r, u
  call check_status(create_material_data_manager(m, b, 2))
  call check_status(material_data_manager_get_state_0(s0, m))
  call check_status(material_data_manager_get_state_1(s1, m))
  call check_status(material_state_manager_set_non_uniform_material_property( &
       s0, 'YoungModulus', E, LOCAL_STORAGE))
  call check_status(material_state_manager_set_uniform_scalar_material_property( &
       s0, 'PoissonRatio', 0d0))
  ! null external state variables of any type
  call check_status(behaviour_get_hypothesis(h, b))
  call check_status(behaviour_get_number_of_external_state_variables(nesv, b))
  do i = 1, nesv
     call check_status(behaviour_get_external_state_variable_name(n, b, i))
     call check_status(behaviour_get_external_state_variable_type(t, b, i))
     call check_status(get_variable_size(vs, h, t))
     call check_status(material_state_manager_set_uniform_external_state_variable( &
          s0, n, zeros(1:vs), LOCAL_STORAGE))
     call check_status(material_state_manager_set_uniform_external_state_variable( &
          s1, n, zeros(1:vs), LOCAL_STORAGE))
  end do
  ! values of v_esv kept by revert
  call check_status(material_state_manager_set_non_uniform_external_state_variable2( &
       s1, 'v_esv', v_esv, LOCAL_STORAGE, NOUPDATE))
  st = material_state_manager_set_non_uniform_external_state_variable( &
       s1, 'v_esv', v_esv(1:3), LOCAL_STORAGE)
  r = check(st % exit_status .ne. MGIS_SUCCESS, 'invalid number of values accepted')
  st = material_state_manager_set_uniform_external_state_variable( &
       s1, 'v_esv', v_esv(1:6:2), LOCAL_STORAGE)
  r = check(st % exit_status .ne. MGIS_SUCCESS, 'non contiguous array accepted')
  call check_status(revert_material_data_manager(m))
  ! integration: v_isv is equal to v_esv and sig_xx to E eto_xx
  call check_status(material_state_manager_get_gradients(g, s1))
  g(1, :) = 1d-3
  call check_status(integrate_material_data_manager_part( &
       ri, m, INTEGRATION_NO_TANGENT_OPERATOR, 1d0, 1, 2))
  r = check(ri == 1, 'integration failed')
  call check_status(material_state_manager_get_thermodynamic_forces(sig, s1))
  call check_status(material_state_manager_get_internal_state_variables(isvs, s1))
  do j = 1, 2
     r = check(abs(sig(1, j) - E(j) * 1d-3) < 1d-8 * E(j) * 1d-3, &
          'invalid stress value')
     do i = 1, 3
        r = check(abs(isvs(i, j) - v_esv(3 * (j - 1) + i)) < 1.d-14, &
             'invalid internal state variable value')
     end do
  end do
  ! mass density per integration point
  call check_status(material_state_manager_set_non_uniform_mass_density( &
       s1, (/ 7800d0, 7900d0 /), LOCAL_STORAGE))
  call check_status(material_state_manager_is_mass_density_uniform(u, s1))
  r = check(.not. u, 'the mass density shall not be uniform')
  ! free ressources
  call check_status(free_material_data_manager(m))
end subroutine check_material_data_manager

program main
  use mgis_behaviour
  use mgis_testing_utilities
  ! behaviour
  type(Behaviour) b
  ! start of the check
  call check_status(load_behaviour(b, &
       get_mfront_behaviour_test_library_path(), &
       'TensorialExternalStateVariableTest', 'Tridimensional'))
  call check_behaviour(b)
  call check_behaviour_data(b)
  call check_material_data_manager(b)
  call tests_summary()
  ! free ressources
  call check_status(free_behaviour(b))
  if(.not. status) then
     stop -1
  end if
end program main

