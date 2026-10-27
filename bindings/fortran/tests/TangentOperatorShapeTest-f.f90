subroutine test()
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  type(Behaviour) b
  type(FiniteStrainBehaviourOptions) o
  type(BehaviourData) d
  type(MaterialDataManager) m
  real(kind=8), dimension (:), pointer :: K1 => null()
  real(kind=8), dimension (:,:), pointer :: K2 => null()
  real(kind=8), dimension (:,:,:), pointer :: K => null()
  integer :: nb
  logical :: r
  call check_status(create_finite_strain_behaviour_options(o))
  call check_status(load_finite_strain_behaviour(b, o, &
       get_mfront_behaviour_test_library_path(), &
       'FiniteStrainSingleCrystal', 'Tridimensional'))
  call check_status(behaviour_get_number_of_tangent_operator_blocks(nb, b))
  r = check(nb == 1, 'invalid number of tangent operator blocks')
  call check_status(allocate_behaviour_data(d, b))
  call check_status(behaviour_data_get_tangent_operator(K2, d))
  r = check((size(K2, 1) == 9) .and. (size(K2, 2) == 6), &
       'invalid shape of the tangent operator (behaviour data)')
  call check_status(behaviour_data_get_tangent_operator_array(K1, d))
  r = check(size(K1) == 54, 'invalid size of the tangent operator array')
  call check_status(create_material_data_manager(m, b, 2))
  call check_status(material_data_manager_allocate_array_of_tangent_operator_blocks(m))
  call check_status(material_data_manager_get_tangent_operator(K, m))
  r = check((size(K, 1) == 9) .and. (size(K, 2) == 6) .and. &
       (size(K, 3) == 2), &
       'invalid shape of the tangent operator (material data manager)')
  call check_status(free_material_data_manager(m))
  call check_status(free_behaviour_data(d))
  call check_status(free_behaviour(b))
  call check_status(free_finite_strain_behaviour_options(o))
end subroutine test

program main
  use mgis_testing_utilities
  call test()
  call tests_summary()
  if(.not. status) then
     stop -1
  end if
end program main
