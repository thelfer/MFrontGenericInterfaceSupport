subroutine test()
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  type(Behaviour) b
  ! size of an array able to hold the gradients' value
  integer gs_size
  ! size of an array able to hold the thermodynamic forces' value
  integer tfs_size
  logical :: r
  ! start of the check
  call check_status(load_behaviour(b, &
       get_mfront_behaviour_test_library_path(), &
       'FiniteStrainSingleCrystal', 'Tridimensional'))
  ! gradients size
  call check_status(behaviour_get_gradients_size(gs_size,b))
  r = check(gs_size==9,'invalid gradients size')
  ! thermodynamic forces size
  call check_status(behaviour_get_thermodynamic_forces_size(tfs_size,b))
  r = check(tfs_size==6,'invalid thermodynamic forces size')
  ! free ressources
  call check_status(free_behaviour(b))
end subroutine test

subroutine test2()
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  type(Behaviour) b
  type(FiniteStrainBehaviourOptions) o
  ! size of an array able to hold the gradients' value
  integer gs_size
  ! size of an array able to hold the thermodynamic forces' value
  integer tfs_size
  logical :: r
  ! start of the check
  call check_status(create_finite_strain_behaviour_options(o))
  call check_status(load_finite_strain_behaviour(b, o, &
       get_mfront_behaviour_test_library_path(), &
       'FiniteStrainSingleCrystal', 'Tridimensional'))
  ! gradients size
  call check_status(behaviour_get_gradients_size(gs_size,b))
  r = check(gs_size==9,'invalid gradients size')
  ! thermodynamic forces size
  call check_status(behaviour_get_thermodynamic_forces_size(tfs_size,b))
  r = check(tfs_size==6,'invalid thermodynamic forces size')
  ! free ressources
  call check_status(free_finite_strain_behaviour_options(o))
  call check_status(free_behaviour(b))
end subroutine test2

subroutine test3(ss, es)
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  character(len=*), intent(in) :: ss
  integer, intent(in) :: es
  type(Behaviour) b
  type(FiniteStrainBehaviourOptions) o
  ! size of an array able to hold the gradients' value
  integer gs_size
  ! size of an array able to hold the thermodynamic forces' value
  integer tfs_size
  logical :: r
  ! start of the check
  call check_status(create_finite_strain_behaviour_options(o))
  call check_status(finite_strain_behaviour_options_set_stress_measure_by_string(o,ss))
  call check_status(load_finite_strain_behaviour(b, o, &
       get_mfront_behaviour_test_library_path(), &
       'FiniteStrainSingleCrystal', 'Tridimensional'))
  ! gradients size
  call check_status(behaviour_get_gradients_size(gs_size,b))
  r = check(gs_size==9,'invalid gradients size')
  ! thermodynamic forces size
  call check_status(behaviour_get_thermodynamic_forces_size(tfs_size,b))
  r = check(tfs_size==es,'invalid thermodynamic forces size')
  ! free ressources
  call check_status(free_finite_strain_behaviour_options(o))
  call check_status(free_behaviour(b))
end subroutine test3

program main
  use mgis_testing_utilities
  call test()
  call test2()
  call test3('PK1',9)
  call test3('PK2',6)
  call test3('CAUCHY',6)
  call tests_summary()
  if(.not. status) then
     stop -1
  end if
end program main
