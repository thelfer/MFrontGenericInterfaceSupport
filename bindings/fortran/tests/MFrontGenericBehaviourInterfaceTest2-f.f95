subroutine check_mp(b, i, e)
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  ! behaviour
  type(behaviour), intent(in) :: b
  ! index of the material property
  integer, intent (in) :: i
  ! expected name of the material property
  character(len=*), intent(in) :: e
  ! name of a material property
  character(len=:), allocatable :: mpn
  logical :: r
  call check_status(behaviour_get_material_property_name(mpn, b, i))
  r = check_string(mpn, e, "invalid material property name '" // mpn // "'")
  call check_status(behaviour_get_material_property_type(mp_t, b, i))
  r = check(mp_t == SCALAR,"invalid type for material property '" // mpn // "'")
end subroutine check_mp

subroutine test()
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  type(Behaviour) b
  logical :: r
  ! library
  character(len=:), allocatable :: l
  ! behaviour name
  character(len=:), allocatable :: bn
  ! behaviour name
  character(len=:), allocatable :: fn
  ! hypothesis
  character(len=:), allocatable :: h
  ! source
  character(len=:), allocatable :: s
  ! version
  character(len=:), allocatable :: v
  ! name of the first internal state variable
  character(len=:), allocatable :: eel
  ! name of the second internal state variable
  character(len=:), allocatable :: p
  ! name of the third internal state variable
  character(len=:), allocatable :: pm
  ! name of the fourth internal state variable
  character(len=:), allocatable :: f
  ! name of the first external state variable
  character(len=:), allocatable :: T
  ! size of an array able to hold the gradients' value
  integer gs_size
  ! size of an array able to hold the thermodynamic forces' value
  integer tfs_size
  ! number of material properties
  integer mps_size
  ! number of material properties
  integer isvs_size
  ! number of external state variables
  integer esvs_size
  ! type of the first internal state variable
  integer eel_t
  ! type of a material property
  integer mp_t
  ! type of the second internal state variable
  integer p_t
  ! type of the third internal state variable
  integer pm_t
  ! type of the fourth internal state variable
  integer f_t
  ! type of the first external state variable
  integer T_t
  ! loop index
  integer i
  ! start of the check
  call check_status(load_behaviour(b, &
       get_mfront_behaviour_test_library_path(), &
       'FiniteStrainSingleCrystal', 'Tridimensional'))
  call check_status(behaviour_get_library(l, b))
  r = check_string(l, get_mfront_behaviour_test_library_path(), &
       'invalid library')
  ! behaviour name
  call check_status(behaviour_get_behaviour_name(bn, b))
  r = check_string(bn, 'FiniteStrainSingleCrystal', 'invalid behaviour name')
  ! function name
  call check_status(behaviour_get_function_name(fn, b))
  r = check_string(fn, 'FiniteStrainSingleCrystal_Tridimensional', &
       'invalid function name')
  ! hypothesis
  call  check_status(behaviour_get_hypothesis(h, b))
  r = check_string(h, 'Tridimensional', 'invalid hypothesis')
  ! source
  call check_status(behaviour_get_source(s, b))
  r = check_string(s, 'FiniteStrainSingleCrystal.mfront', 'invalid source')
  ! version
  call check_status(behaviour_get_tfel_version(v, b))
  r = check_string(v, get_tfel_version(), "invalid TFEL version")
  ! gradients size
  call check_status(behaviour_get_gradients_size(gs_size,b))
  r = check(gs_size==9,'invalid gradients size')
  ! thermodynamic forces size
  call check_status(behaviour_get_thermodynamic_forces_size(tfs_size,b))
  r = check(tfs_size==6,'invalid thermodynamic forces size')
  ! material properties
  call check_status(behaviour_get_number_of_material_properties(mps_size, b))
  if (check(mps_size .eq. 16, "invalid number of material properties")) then
     call check_mp(b,1,'YoungModulus1')
     call check_mp(b,2,'YoungModulus2')
     call check_mp(b,3,'YoungModulus3')
     call check_mp(b,4,'PoissonRatio12')
     call check_mp(b,5,'PoissonRatio23')
     call check_mp(b,6,'PoissonRatio13')
     call check_mp(b,7,'ShearModulus12')
     call check_mp(b,8,'ShearModulus23')
     call check_mp(b,9,'ShearModulus13')
     call check_mp(b,10,'m')
     call check_mp(b,11,'K')
     call check_mp(b,12,'C')
     call check_mp(b,13,'R0')
     call check_mp(b,14,'Q')
     call check_mp(b,15,'b')
     call check_mp(b,16,'d1')
  end if
  ! ! internal state variable
  ! call check_status(behaviour_get_number_of_internal_state_variables(isvs_size, b))
  ! if (check(isvs_size .eq. 4, "invalid number of internal state variable")) then
  !    call check_status(behaviour_get_internal_state_variable_name(eel, b, 0))
  !    r = check_string(eel, "ElasticStrain", "invalid internal state variable name")
  !    call check_status(behaviour_get_internal_state_variable_type(eel_t, b, 0))
  !    r = check(eel_t == STENSOR,"invalid type for internal state variable 'ElasticStrain'")
  !    call check_status(behaviour_get_internal_state_variable_name(p, b, 1))
  !    r = check_string(p, "EquivalentPlasticStrain", "invalid internal state variable name")
  !    call check_status(behaviour_get_internal_state_variable_type(p_t, b, 1))
  !    r = check(p_t == SCALAR,"invalid type for internal state variable 'EquivalentPlasticStrain'")
  !    call check_status(behaviour_get_internal_state_variable_name(pm, b, 2))
  !    r = check_string(pm, "MatrixEquivalentPlasticStrain", "invalid internal state variable name")
  !    call check_status(behaviour_get_internal_state_variable_type(pm_t, b, 2))
  !    r = check(pm_t == SCALAR,"invalid type for internal state variable 'MatrixEquivalentPlasticStrain'")
  !    call check_status(behaviour_get_internal_state_variable_name(f, b, 3))
  !    r = check_string(f, "Porosity", "invalid internal state variable name")
  !    call check_status(behaviour_get_internal_state_variable_type(f_t, b, 3))
  !    r = check(f_t == SCALAR,"invalid type for internal state variable 'Porosity'")
  ! end if
  ! external state variable
  call check_status(behaviour_get_number_of_external_state_variables(esvs_size, b))
  if (check(esvs_size .eq. 1, "invalid number of external state variable")) then
     call check_status(behaviour_get_external_state_variable_name(T, b, 1))
     r = check_string(T, "Temperature", "invalid external state variable name")
     call check_status(behaviour_get_external_state_variable_type(T_t, b, 1))
     r = check(T_t == SCALAR,"invalid type for external state variable 'Temperature'")
  end if
  ! free ressources
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
