subroutine test_initialize_function()
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  character(len=*), parameter :: f = 'StressFromInitialPressure'
  real(kind=8), parameter :: pr = -1.2d5
  real(kind=8), parameter :: eps = -10 * pr * epsilon(pr)
  type(ThreadPool) p
  type(Behaviour) b
  type(MaterialDataManager) m
  type(MaterialStateManager) s0, s1
  real(kind=8), dimension(:,:), pointer :: sig => null()
  character(len=:), allocatable :: n
  integer c, t, ri, i, j
  logical :: r
  call check_status(load_behaviour(b, &
       get_mfront_behaviour_test_library_path(), &
       'InitializeFunctionTest', 'Tridimensional'))
  ! description
  call check_status(behaviour_get_number_of_initialize_functions(c, b))
  r = check(c == 3, 'invalid number of initialize functions')
  call check_status(behaviour_get_initialize_function_name(n, b, 3))
  r = check_string(n, f, 'invalid initialize function name')
  call check_status(behaviour_get_initialize_function_inputs_size(c, b, f))
  r = check(c == 1, 'invalid size of the inputs')
  call check_status(behaviour_get_number_of_initialize_function_inputs(c, b, f))
  r = check(c == 1, 'invalid number of inputs')
  call check_status(behaviour_get_initialize_function_input_name(n, b, f, 1))
  r = check_string(n, 'pr', 'invalid input name')
  call check_status(behaviour_get_initialize_function_input_type(t, b, f, 1))
  r = check(t == SCALAR, 'invalid input type')
  call check_status(behaviour_get_initialize_function_input_offset(c, b, f, 1))
  r = check(c == 1, 'invalid input offset')
  ! execution on a material data manager, using a thread pool
  call check_status(create_thread_pool(p, 2))
  call check_status(create_material_data_manager(m, b, 2))
  call check_status(material_data_manager_get_state_0(s0, m))
  call check_status(material_data_manager_get_state_1(s1, m))
  call check_status(material_state_manager_set_us_external_state_variable( &
       s0, 'Temperature', 293.15d0))
  call check_status(material_state_manager_set_us_external_state_variable( &
       s1, 'Temperature', 293.15d0))
  call check_status(execute_initialize_function_material_data_manager( &
       ri, p, m, f, [pr]))
  r = check(ri == 1, 'invalid result of the initialize function')
  call check_status(material_state_manager_get_thermodynamic_forces(sig, s1))
  do j = 1, 2
     do i = 1, 6
        if (i <= 3) then
           r = check(abs(sig(i, j) - pr) < eps, 'invalid stress value')
        else
           r = check(abs(sig(i, j)) < eps, 'invalid stress value')
        end if
     end do
  end do
  ! free ressources
  call check_status(free_material_data_manager(m))
  call check_status(free_behaviour(b))
  call check_status(free_thread_pool(p))
end subroutine test_initialize_function

subroutine test_post_processing()
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  character(len=*), parameter :: pp = 'PrincipalStrain'
  real(kind=8), dimension(3), parameter :: e2 = [1.2d-2, 1.3d-2, 1.4d-2]
  real(kind=8), parameter :: eps = 10 * epsilon(1.d0)
  type(Behaviour) b
  type(BehaviourData) d
  type(State) s0, s1
  real(kind=8), dimension(6) :: e
  real(kind=8), dimension(3) :: outputs
  character(len=:), allocatable :: n
  integer c, t, ri, i
  logical :: r
  call check_status(load_behaviour(b, &
       get_mfront_behaviour_test_library_path(), &
       'PostProcessingTest', 'Tridimensional'))
  ! description
  call check_status(behaviour_get_number_of_post_processings(c, b))
  r = check(c == 1, 'invalid number of post-processings')
  call check_status(behaviour_get_post_processing_name(n, b, 1))
  r = check_string(n, pp, 'invalid post-processing name')
  call check_status(behaviour_get_post_processing_outputs_size(c, b, pp))
  r = check(c == 3, 'invalid size of the outputs')
  call check_status(behaviour_get_number_of_post_processing_outputs(c, b, pp))
  r = check(c == 1, 'invalid number of outputs')
  call check_status(behaviour_get_post_processing_output_name(n, b, pp, 1))
  r = check_string(n, 'PrincipalStrain', 'invalid output name')
  call check_status(behaviour_get_post_processing_output_type(t, b, pp, 1))
  r = check(t == VECTOR_3D, 'invalid output type')
  call check_status(behaviour_get_post_processing_output_offset(c, b, pp, 1))
  r = check(c == 1, 'invalid output offset')
  ! execution on a behaviour data
  call check_status(allocate_behaviour_data(d, b))
  call check_status(behaviour_data_get_state_0(s0, d))
  call check_status(behaviour_data_get_state_1(s1, d))
  call check_status(state_set_scalar_external_state_variable_by_name( &
       s0, 'Temperature', 293.15d0))
  call check_status(state_set_scalar_external_state_variable_by_name( &
       s1, 'Temperature', 293.15d0))
  e = (/ 1.3d-2, 1.2d-2, 1.4d-2, 0d0, 0d0, 0d0 /)
  call check_status(state_set_gradient_by_name(s1, 'Strain', e))
  call check_status(execute_post_processing(ri, outputs, d, b, pp))
  r = check(ri == 0, 'invalid result of the post-processing')
  do i = 1, 3
     r = check(abs(outputs(i) - e2(i)) < eps, 'invalid output value')
  end do
  ! free ressources
  call check_status(free_behaviour_data(d))
  call check_status(free_behaviour(b))
end subroutine test_post_processing

program main
  use mgis_testing_utilities
  call test_initialize_function()
  call test_post_processing()
  call tests_summary()
  if(.not. status) then
     stop -1
  end if
end program main
