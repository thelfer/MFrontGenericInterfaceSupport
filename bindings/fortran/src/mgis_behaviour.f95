module mgis_behaviour
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
  enum, bind(C)
     enumerator :: SCALAR  = 0
     enumerator :: VECTOR  = 1
     enumerator :: STENSOR = 2
     enumerator :: TENSOR  = 3
  end enum
  enum, bind(C)
     enumerator :: ISOTROPIC   = 0
     enumerator :: ORTHOTROPIC = 1
  end enum
  enum, bind(C)
     enumerator :: GENERALBEHAVIOUR = 0
     enumerator :: STANDARDSTRAINBASEDBEHAVIOUR = 1
     enumerator :: STANDARDFINITESTRAINBEHAVIOUR = 2
     enumerator :: COHESIVEZONEMODEL = 3
  end enum
  enum, bind(C)
     enumerator :: UNDEFINEDKINEMATIC = 0
     enumerator :: SMALLSTRAINKINEMATIC = 1
     enumerator :: COHESIVEZONEKINEMATIC = 2
     enumerator :: FINITESTRAINKINEMATIC_F_CAUCHY = 3
     enumerator :: FINITESTRAINKINEMATIC_ETO_PK1 = 4
  end enum
  enum, bind(C)
     enumerator :: LOCAL_STORAGE    = 0
     enumerator :: EXTERNAL_STORAGE = 1
  end enum
  enum, bind(C)
     enumerator :: PREDICTION_TANGENT_OPERATOR = -3
     enumerator :: PREDICTION_SECANT_OPERATOR = -2
     enumerator :: PREDICTION_ELASTIC_OPERATOR = -1
     enumerator :: INTEGRATION_NO_TANGENT_OPERATOR = 0
     enumerator :: INTEGRATION_ELASTIC_OPERATOR = 1
     enumerator :: INTEGRATION_SECANT_OPERATOR = 2
     enumerator :: INTEGRATION_TANGENT_OPERATOR = 3
     enumerator :: INTEGRATION_CONSISTENT_TANGENT_OPERATOR = 4
  end enum
  enum, bind(C)
     enumerator :: CAUCHY = 0
     enumerator :: PK2 = 1
     enumerator :: PK1 = 2
  end enum
  enum, bind(C)
     enumerator :: DSIG_DF = 0
     enumerator :: DS_DEGL = 1
     enumerator :: DPK1_DF = 2
  end enum
  type :: FiniteStrainBehaviourOptions
    private
    type(c_ptr) :: ptr = c_null_ptr
  end type FiniteStrainBehaviourOptions
  type :: Behaviour
    private
    type(c_ptr) :: ptr = c_null_ptr
  end type Behaviour
  type :: BehaviourData
    private
    type(c_ptr) :: ptr = c_null_ptr
  end type BehaviourData
  type :: State
    private
    type(c_ptr) :: ptr = c_null_ptr
  end type State
  type :: MaterialStateManagerInitializer
    private
    type(c_ptr) :: ptr = c_null_ptr
  end type MaterialStateManagerInitializer
  type :: MaterialStateManager
    private
    type(c_ptr) :: ptr = c_null_ptr
  end type MaterialStateManager
  type :: MaterialDataManagerInitializer
    private
    type(c_ptr) :: ptr = c_null_ptr
  end type MaterialDataManagerInitializer
  type :: MaterialDataManager
    private
    type(c_ptr) :: ptr = c_null_ptr
  end type MaterialDataManager
contains
  !
  function get_space_dimension(vs, h) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t, c_int
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function get_space_dimension_wrapper(vs,h) &
            bind(c,name = 'mgis_bv_get_space_dimension') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_size_t, c_char
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_size_t), intent(out) :: vs
         character(len=1,kind=c_char), dimension(*), intent(in) :: h
         type(mgis_status) :: r
       end function get_space_dimension_wrapper
    end interface
    integer(kind=c_int), intent(out) :: vs
    character(len=*), intent(in) :: h
    type(mgis_status) :: s
    integer(kind=c_size_t) :: ns
    s = get_space_dimension_wrapper(ns, convert_fortran_string(h))
    vs = ns
  end function get_space_dimension
  !
  function get_stensor_size(vs, h) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t, c_int
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function get_stensor_size_wrapper(vs,h) &
            bind(c,name = 'mgis_bv_get_stensor_size') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_size_t, c_char
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_size_t), intent(out) :: vs
         character(len=1,kind=c_char), dimension(*), intent(in) :: h
         type(mgis_status) :: r
       end function get_stensor_size_wrapper
    end interface
    integer(kind=c_int), intent(out) :: vs
    character(len=*), intent(in) :: h
    type(mgis_status) :: s
    integer(kind=c_size_t) :: ns
    s = get_stensor_size_wrapper(ns, convert_fortran_string(h))
    vs = ns
  end function get_stensor_size
  !
  function get_tensor_size(vs, h) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t, c_int
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function get_tensor_size_wrapper(vs,h) &
            bind(c,name = 'mgis_bv_get_tensor_size') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_size_t, c_char
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_size_t), intent(out) :: vs
         character(len=1,kind=c_char), dimension(*), intent(in) :: h
         type(mgis_status) :: r
       end function get_tensor_size_wrapper
    end interface
    integer(kind=c_int), intent(out) :: vs
    character(len=*), intent(in) :: h
    type(mgis_status) :: s
    integer(kind=c_size_t) :: ns
    s = get_tensor_size_wrapper(ns, convert_fortran_string(h))
    vs = ns
  end function get_tensor_size
  !
  function get_variable_size(vs, h, t) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t, c_int
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function get_variable_size_wrapper(vs,h,t) &
            bind(c,name = 'mgis_bv_get_variable_size') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_size_t, c_char, c_int
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_size_t), intent(out) :: vs
         character(len=1,kind=c_char), dimension(*), intent(in) :: h
         integer(kind=c_int), intent(in), value :: t
         type(mgis_status) :: r
       end function get_variable_size_wrapper
    end interface
    integer(kind=c_int), intent(out) :: vs
    character(len=*), intent(in) :: h
    integer(kind=c_int), intent(in) :: t
    type(mgis_status) :: s
    integer(kind=c_size_t) :: ns
    s = get_variable_size_wrapper(ns, convert_fortran_string(h), t)
    vs = ns
  end function get_variable_size
  !
  function create_finite_strain_behaviour_options(o) result(s)
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function create_finite_strain_behaviour_options_wrapper(ptr) &
            bind(c,name = 'mgis_bv_create_finite_strain_behaviour_options') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: ptr
         type(mgis_status) :: r
       end function create_finite_strain_behaviour_options_wrapper
    end interface
    type(FiniteStrainBehaviourOptions), intent(out) :: o
    type(mgis_status) :: s
    s = create_finite_strain_behaviour_options_wrapper(o%ptr)
  end function create_finite_strain_behaviour_options
  !
  function finite_strain_behaviour_options_set_stress_measure(o, ss) result(s)
    use, intrinsic :: iso_c_binding, only: c_int
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function fsb_opts_set_stress_measure_wrapper(o,ss) &
            bind(c,name = 'mgis_bv_finite_strain_behaviour_options_set_stress_measure') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_int
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(in), value :: o
         integer(kind=c_int), intent(in),value :: ss
         type(mgis_status) :: r
       end function fsb_opts_set_stress_measure_wrapper
    end interface
    type(FiniteStrainBehaviourOptions), intent(in) :: o
    integer(kind=c_int), intent(in) :: ss
    type(mgis_status) :: s
    s = fsb_opts_set_stress_measure_wrapper(o%ptr,ss)
  end function finite_strain_behaviour_options_set_stress_measure
  !
  function finite_strain_behaviour_options_set_stress_measure_by_string(o, ss) result(s)
    use, intrinsic :: iso_c_binding, only: c_int
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function fsb_opts_set_stress_measure_by_string_wrapper(o,ss) &
            bind(c,name = 'mgis_bv_finite_strain_behaviour_options_set_stress_measure_by_string') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_char
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(in), value :: o
         character(len=1,kind=c_char), dimension(*), intent(in) :: ss
         type(mgis_status) :: r
       end function fsb_opts_set_stress_measure_by_string_wrapper
    end interface
    type(FiniteStrainBehaviourOptions), intent(in) :: o
    character(len=*), intent(in) :: ss
    type(mgis_status) :: s
    s = fsb_opts_set_stress_measure_by_string_wrapper(o%ptr, convert_fortran_string(ss))
  end function finite_strain_behaviour_options_set_stress_measure_by_string
  !
  function finite_strain_behaviour_options_set_tangent_operator(o, to) result(s)
    use, intrinsic :: iso_c_binding, only: c_int
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function fsb_opts_set_tangent_operator_wrapper(o,to) &
            bind(c,name = 'mgis_bv_finite_strain_behaviour_options_set_tangent_operator') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_int
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(in), value :: o
         integer(kind=c_int), intent(in),value :: to
         type(mgis_status) :: r
       end function fsb_opts_set_tangent_operator_wrapper
    end interface
    type(FiniteStrainBehaviourOptions), intent(in) :: o
    integer(kind=c_int), intent(in) :: to
    type(mgis_status) :: s
    s = fsb_opts_set_tangent_operator_wrapper(o%ptr,to)
  end function finite_strain_behaviour_options_set_tangent_operator
  !
  function finite_strain_behaviour_options_set_tangent_operator_by_string(o, to) result(s)
    use, intrinsic :: iso_c_binding, only: c_int
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function fsb_opts_set_tangent_operator_by_string_wrapper(o,to) &
            bind(c,name = 'mgis_bv_finite_strain_behaviour_options_set_tangent_operator_by_string') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_char
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(in), value :: o
         character(len=1,kind=c_char), dimension(*), intent(in) :: to
         type(mgis_status) :: r
       end function fsb_opts_set_tangent_operator_by_string_wrapper
    end interface
    type(FiniteStrainBehaviourOptions), intent(in) :: o
    character(len=*), intent(in) :: to
    type(mgis_status) :: s
    s = fsb_opts_set_tangent_operator_by_string_wrapper(o%ptr,convert_fortran_string(to))
  end function finite_strain_behaviour_options_set_tangent_operator_by_string
  !
  function free_finite_strain_behaviour_options(ptr) result(r)
    use, intrinsic :: iso_c_binding, only: c_associated
    use mgis
    implicit none
    interface
       function free_finite_strain_behaviour_options_wrapper(ptr) &
            bind(c, name='mgis_bv_free_finite_strain_behaviour_options') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis
         implicit none
         type(c_ptr), intent(inout) :: ptr
         type(mgis_status) :: r
       end function free_finite_strain_behaviour_options_wrapper
    end interface
    type(FiniteStrainBehaviourOptions), intent(inout) :: ptr
    type(mgis_status) :: r
    if (c_associated(ptr%ptr)) then
       r = free_finite_strain_behaviour_options_wrapper(ptr%ptr)
    end if
  end function free_finite_strain_behaviour_options
  !
  function is_standard_finite_strain_behaviour(b,l,bn) result(s)
    use, intrinsic :: iso_c_binding, only: c_int
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    implicit none
    interface
       function is_standard_finite_strain_behaviour_wrapper(b,l,bn) &
            bind(c,name = 'mgis_bv_is_standard_finite_strain_behaviour') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_char, c_ptr, c_int
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_int), intent(out) :: b
         character(len=1,kind=c_char), dimension(*), intent(in) :: l
         character(len=1,kind=c_char), dimension(*), intent(in) :: bn
         type(mgis_status) :: r
       end function is_standard_finite_strain_behaviour_wrapper
    end interface
    logical, intent(out) :: b
    character(len=*), intent(in) :: l
    character(len=*), intent(in) :: bn
    type(mgis_status) :: s
    integer(kind=c_int) :: r
    s = is_standard_finite_strain_behaviour_wrapper(r, convert_fortran_string(l), &
         convert_fortran_string(bn))
    if (r .eq. 0) then
       b = .false.
    else
       b = .true.
    end if
  end function is_standard_finite_strain_behaviour
  !
  function load_behaviour(b,l,bn,h) result(s)
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    implicit none
    interface
       function load_behaviour_wrapper(ptr,l,b,h) &
            bind(c,name = 'mgis_bv_load_behaviour') result(r)
         use, intrinsic :: iso_c_binding, only: c_char, c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: ptr
         character(len=1,kind=c_char), dimension(*), intent(in) :: l
         character(len=1,kind=c_char), dimension(*), intent(in) :: b
         character(len=1,kind=c_char), dimension(*), intent(in) :: h
         type(mgis_status) :: r
       end function load_behaviour_wrapper
    end interface
    type(behaviour), intent(out) :: b
    character(len=*), intent(in) :: l
    character(len=*), intent(in) :: bn
    character(len=*), intent(in) :: h
    type(mgis_status) :: s
    s = load_behaviour_wrapper(b%ptr, convert_fortran_string(l), &
         convert_fortran_string(bn), &
         convert_fortran_string(h))
  end function load_behaviour
  !
  function load_finite_strain_behaviour(b,o,l,bn,h) result(s)
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    implicit none
    interface
       function load_finite_strain_behaviour_wrapper(ptr,o,l,b,h) &
            bind(c,name = 'mgis_bv_load_finite_strain_behaviour') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_char, c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: ptr
         type(c_ptr), intent(in), value :: o
         character(len=1,kind=c_char), dimension(*), intent(in) :: l
         character(len=1,kind=c_char), dimension(*), intent(in) :: b
         character(len=1,kind=c_char), dimension(*), intent(in) :: h
         type(mgis_status) :: r
       end function load_finite_strain_behaviour_wrapper
    end interface
    type(Behaviour), intent(out) :: b
    type(FiniteStrainBehaviourOptions), intent(in) :: o
    character(len=*), intent(in) :: l
    character(len=*), intent(in) :: bn
    character(len=*), intent(in) :: h
    type(mgis_status) :: s
    s = load_finite_strain_behaviour_wrapper(b%ptr, o%ptr, &
         convert_fortran_string(l), &
         convert_fortran_string(bn), &
         convert_fortran_string(h))
  end function load_finite_strain_behaviour
  ! behaviour_get_tangent_operator_array_size
  function behaviour_get_tangent_operator_array_size(n,b) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    implicit none
    interface
       function behaviour_get_tangent_operator_array_size_wrapper(l,b) &
            bind(c,name = 'mgis_bv_behaviour_get_tangent_operator_array_size') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_size_t, c_ptr
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_size_t), intent(out) :: l
         type(c_ptr), intent(in), value :: b
         type(mgis_status) :: r
       end function behaviour_get_tangent_operator_array_size_wrapper
    end interface
    integer :: n
    type(behaviour), intent(in) :: b
    type(mgis_status) :: s
    integer(kind=c_size_t) :: ns
    s = behaviour_get_tangent_operator_array_size_wrapper(ns, b%ptr)
    n = ns
  end function behaviour_get_tangent_operator_array_size
  ! behaviour_get_library
  function behaviour_get_library(l,b) result(s)
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_get_library_wrapper(l,b) bind(c,name = 'mgis_bv_behaviour_get_library') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: l
         type(c_ptr), intent(in), value :: b
         type(mgis_status) :: r
       end function behaviour_get_library_wrapper
    end interface
    character(len=:), allocatable, intent(out) :: l
    type(behaviour), intent(in) :: b
    type(mgis_status) :: s
    type(c_ptr) :: lp
    s = behaviour_get_library_wrapper(lp,b%ptr)
    if (s % exit_status == MGIS_SUCCESS) then
       l = convert_c_string(lp)
    else
       l = get_empty_string();
    end if
  end function behaviour_get_library
  ! behaviour_get_behaviour_name
  function behaviour_get_behaviour_name(l,b) result(s)
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_get_behaviour_name_wrapper(l,b) bind(c,name = 'mgis_bv_behaviour_get_behaviour_name') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: l
         type(c_ptr), intent(in), value :: b
         type(mgis_status) :: r
       end function behaviour_get_behaviour_name_wrapper
    end interface
    character(len=:), allocatable, intent(out) :: l
    type(behaviour), intent(in) :: b
    type(mgis_status) :: s
    type(c_ptr) :: lp
    s = behaviour_get_behaviour_name_wrapper(lp, b%ptr)
    if (s % exit_status == MGIS_SUCCESS) then
       l = convert_c_string(lp)
    else
       l = get_empty_string();
    end if
  end function behaviour_get_behaviour_name
  ! behaviour_get_function_name
  function behaviour_get_function_name(l,b) result(s)
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_get_function_name_wrapper(l,b) bind(c,name = 'mgis_bv_behaviour_get_function_name') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: l
         type(c_ptr), intent(in), value :: b
         type(mgis_status) :: r
       end function behaviour_get_function_name_wrapper
    end interface
    character(len=:), allocatable, intent(out) :: l
    type(behaviour), intent(in) :: b
    type(mgis_status) :: s
    type(c_ptr) :: lp
    s = behaviour_get_function_name_wrapper(lp, b%ptr)
    if (s % exit_status == MGIS_SUCCESS) then
       l = convert_c_string(lp)
    else
       l = get_empty_string();
    end if
  end function behaviour_get_function_name
  ! behaviour_get_source
  function behaviour_get_source(l,b) result(s)
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_get_source_wrapper(l,b) bind(c,name = 'mgis_bv_behaviour_get_source') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: l
         type(c_ptr), intent(in), value :: b
         type(mgis_status) :: r
       end function behaviour_get_source_wrapper
    end interface
    character(len=:), allocatable, intent(out) :: l
    type(behaviour), intent(in) :: b
    type(mgis_status) :: s
    type(c_ptr) :: lp
    s = behaviour_get_source_wrapper(lp, b%ptr)
    if (s % exit_status == MGIS_SUCCESS) then
       l = convert_c_string(lp)
    else
       l = get_empty_string();
    end if
  end function behaviour_get_source
  ! behaviour_get_hypothesis
  function behaviour_get_hypothesis(l,b) result(s)
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_get_hypothesis_wrapper(l,b) bind(c,name = 'mgis_bv_behaviour_get_hypothesis') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: l
         type(c_ptr), intent(in), value :: b
         type(mgis_status) :: r
       end function behaviour_get_hypothesis_wrapper
    end interface
    character(len=:), allocatable, intent(out) :: l
    type(behaviour), intent(in) :: b
    type(mgis_status) :: s
    type(c_ptr) :: lp
    s = behaviour_get_hypothesis_wrapper(lp, b%ptr)
    if (s % exit_status == MGIS_SUCCESS) then
       l = convert_c_string(lp)
    else
       l = get_empty_string();
    end if
  end function behaviour_get_hypothesis
  ! behaviour_get_tfel_version
  function behaviour_get_tfel_version(l,b) result(s)
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_get_tfel_version_wrapper(l,b) bind(c,name = 'mgis_bv_behaviour_get_tfel_version') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: l
         type(c_ptr), intent(in), value :: b
         type(mgis_status) :: r
       end function behaviour_get_tfel_version_wrapper
    end interface
    character(len=:), allocatable, intent(out) :: l
    type(behaviour), intent(in) :: b
    type(mgis_status) :: s
    type(c_ptr) :: lp
    s = behaviour_get_tfel_version_wrapper(lp, b%ptr)
    if (s % exit_status == MGIS_SUCCESS) then
       l = convert_c_string(lp)
    else
       l = get_empty_string();
    end if
  end function behaviour_get_tfel_version
  !
  function behaviour_get_gradients_size(n,b) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    implicit none
    interface
       function behaviour_get_gradients_size_wrapper(l,b) &
            bind(c,name = 'mgis_bv_behaviour_get_gradients_size') result(r)
         use, intrinsic :: iso_c_binding, only: c_size_t, c_ptr
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_size_t), intent(out) :: l
         type(c_ptr), intent(in), value :: b
         type(mgis_status) :: r
       end function behaviour_get_gradients_size_wrapper
    end interface
    integer :: n
    type(behaviour), intent(in) :: b
    type(mgis_status) :: s
    integer(kind=c_size_t) :: ns
    s = behaviour_get_gradients_size_wrapper(ns, b%ptr)
    n = ns
  end function behaviour_get_gradients_size
  !
  function behaviour_get_thermodynamic_forces_size(n,b) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    implicit none
    interface
       function behaviour_get_thermodynamic_forces_size_wrapper(l,b) &
            bind(c,name = 'mgis_bv_behaviour_get_thermodynamic_forces_size') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_size_t, c_ptr
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_size_t), intent(out) :: l
         type(c_ptr), intent(in), value :: b
         type(mgis_status) :: r
       end function behaviour_get_thermodynamic_forces_size_wrapper
    end interface
    integer :: n
    type(behaviour), intent(in) :: b
    type(mgis_status) :: s
    integer(kind=c_size_t) :: ns
    s = behaviour_get_thermodynamic_forces_size_wrapper(ns, b%ptr)
    n = ns
  end function behaviour_get_thermodynamic_forces_size
  !
  function behaviour_get_number_of_parameters(n,b) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    implicit none
    interface
       function behaviour_get_number_of_parameters_wrapper(l,b) &
            bind(c,name = 'mgis_bv_behaviour_get_number_of_parameters') result(r)
         use, intrinsic :: iso_c_binding, only: c_size_t, c_ptr
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_size_t), intent(out) :: l
         type(c_ptr), intent(in), value :: b
         type(mgis_status) :: r
       end function behaviour_get_number_of_parameters_wrapper
    end interface
    integer :: n
    type(behaviour), intent(in) :: b
    type(mgis_status) :: s
    integer(kind=c_size_t) :: ns
    s = behaviour_get_number_of_parameters_wrapper(ns, b%ptr)
    n = ns
  end function behaviour_get_number_of_parameters
  !
  function behaviour_get_parameter_name(l, b, n) result(s)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_get_parameter_name_wrapper(l, b, n) &
            bind(c,name = 'mgis_bv_behaviour_get_parameter_name') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: l
         type(c_ptr), intent(in), value :: b
         integer(kind=c_size_t), intent(in), value :: n
         type(mgis_status) :: r
       end function behaviour_get_parameter_name_wrapper
    end interface
    character(len=:), allocatable, intent(out) :: l
    type(behaviour), intent(in) :: b
    integer, intent (in) :: n
    type(mgis_status) :: s
    type(c_ptr) :: lp
    integer(kind = c_size_t) :: nc
    if(.not. convert_to_c_index(nc, n)) then
       s = report_failure("invalid index")
       return
    end if
    s = behaviour_get_parameter_name_wrapper(lp, b%ptr, nc)
    if (s % exit_status == MGIS_SUCCESS) then
       l = convert_c_string(lp)
    else
       l = get_empty_string();
    end if
  end function behaviour_get_parameter_name
  !
  function behaviour_get_parameter_default_value(v, b, n) result(r)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_get_parameter_default_value_wrapper(v, b, n) &
            bind(c,name = 'mgis_bv_behaviour_get_parameter_default_value') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_double
         use mgis, only: mgis_status
         implicit none
         real(kind=c_double), intent(out) :: v
         type(c_ptr), intent(in), value :: b
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: r
       end function behaviour_get_parameter_default_value_wrapper
    end interface
    real(kind=8), intent(out) :: v
    type(behaviour), intent(in) :: b
    character(len=*), intent(in) :: n
    type(mgis_status) :: r
    r = behaviour_get_parameter_default_value_wrapper(&
         v, b%ptr, convert_fortran_string(n))
  end function behaviour_get_parameter_default_value
  !
  function behaviour_get_number_of_material_properties(n,b) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    implicit none
    interface
       function behaviour_get_number_of_material_properties_wrapper(l,b) &
            bind(c,name = 'mgis_bv_behaviour_get_number_of_material_properties') result(r)
         use, intrinsic :: iso_c_binding, only: c_size_t, c_ptr
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_size_t), intent(out) :: l
         type(c_ptr), intent(in), value :: b
         type(mgis_status) :: r
       end function behaviour_get_number_of_material_properties_wrapper
    end interface
    integer :: n
    type(behaviour), intent(in) :: b
    type(mgis_status) :: s
    integer(kind=c_size_t) :: ns
    s = behaviour_get_number_of_material_properties_wrapper(ns, b%ptr)
    n = ns
  end function behaviour_get_number_of_material_properties
  !
  function behaviour_get_material_property_name(l, b, n) result(s)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_get_material_property_name_wrapper(l, b, n) &
            bind(c,name = 'mgis_bv_behaviour_get_material_property_name') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: l
         type(c_ptr), intent(in), value :: b
         integer(kind=c_size_t), intent(in), value :: n
         type(mgis_status) :: r
       end function behaviour_get_material_property_name_wrapper
    end interface
    character(len=:), allocatable, intent(out) :: l
    type(behaviour), intent(in) :: b
    integer, intent (in) :: n
    type(mgis_status) :: s
    type(c_ptr) :: lp
    integer(kind = c_size_t) :: nc
    if(.not. convert_to_c_index(nc, n)) then
       s = report_failure("invalid index")
       return
    end if
    s = behaviour_get_material_property_name_wrapper(lp, b%ptr, nc)
    if (s % exit_status == MGIS_SUCCESS) then
       l = convert_c_string(lp)
    else
       l = get_empty_string();
    end if
  end function behaviour_get_material_property_name
  !
  function behaviour_get_material_property_type(t, b, n) result(s)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function behaviour_get_material_property_type_wrapper(t, b, n) &
            bind(c,name = 'mgis_bv_behaviour_get_material_property_type') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         implicit none
         integer t
         type(c_ptr), intent(in), value :: b
         integer(kind=c_size_t), intent(in), value :: n
         type(mgis_status) :: r
       end function behaviour_get_material_property_type_wrapper
    end interface
    integer, intent (out) :: t
    type(behaviour), intent(in) :: b
    integer, intent (in) :: n
    type(mgis_status) :: s
    integer(kind = c_size_t) :: nc
    if(.not. convert_to_c_index(nc, n)) then
       s = report_failure("invalid index")
       return
    end if
    s = behaviour_get_material_property_type_wrapper(t, b%ptr, nc)
  end function behaviour_get_material_property_type
  !
  function behaviour_get_internal_state_variables_size(n,b) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    implicit none
    interface
       function behaviour_get_internal_state_variables_size_wrapper(l,b) &
            bind(c,name = 'mgis_bv_behaviour_get_internal_state_variables_size') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_size_t, c_ptr
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_size_t), intent(out) :: l
         type(c_ptr), intent(in), value :: b
         type(mgis_status) :: r
       end function behaviour_get_internal_state_variables_size_wrapper
    end interface
    integer :: n
    type(behaviour), intent(in) :: b
    type(mgis_status) :: s
    integer(kind=c_size_t) :: ns
    s = behaviour_get_internal_state_variables_size_wrapper(ns, b%ptr)
    n = ns
  end function behaviour_get_internal_state_variables_size
  !
  function behaviour_get_number_of_internal_state_variables(n,b) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    implicit none
    interface
       function behaviour_get_number_of_internal_state_variables_wrapper(l,b) &
            bind(c,name = 'mgis_bv_behaviour_get_number_of_internal_state_variables') result(r)
         use, intrinsic :: iso_c_binding, only: c_size_t, c_ptr
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_size_t), intent(out) :: l
         type(c_ptr), intent(in), value :: b
         type(mgis_status) :: r
       end function behaviour_get_number_of_internal_state_variables_wrapper
    end interface
    integer :: n
    type(behaviour), intent(in) :: b
    type(mgis_status) :: s
    integer(kind=c_size_t) :: ns
    s = behaviour_get_number_of_internal_state_variables_wrapper(ns, b%ptr)
    n = ns
  end function behaviour_get_number_of_internal_state_variables
  !
  function behaviour_get_internal_state_variable_name(l, b, n) result(s)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_get_internal_state_variable_name_wrapper(l, b, n) &
            bind(c,name = 'mgis_bv_behaviour_get_internal_state_variable_name') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: l
         type(c_ptr), intent(in), value :: b
         integer(kind=c_size_t), intent(in), value :: n
         type(mgis_status) :: r
       end function behaviour_get_internal_state_variable_name_wrapper
    end interface
    character(len=:), allocatable, intent(out) :: l
    type(behaviour), intent(in) :: b
    integer, intent (in) :: n
    type(mgis_status) :: s
    type(c_ptr) :: lp
    integer(kind = c_size_t) :: nc
    if(.not. convert_to_c_index(nc, n)) then
       s = report_failure("invalid index")
       return
    end if
    s = behaviour_get_internal_state_variable_name_wrapper(lp, b%ptr, nc)
    if (s % exit_status == MGIS_SUCCESS) then
       l = convert_c_string(lp)
    else
       l = get_empty_string();
    end if
  end function behaviour_get_internal_state_variable_name
  !
  function behaviour_get_internal_state_variable_type(t, b, n) result(s)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function behaviour_get_internal_state_variable_type_wrapper(t, b, n) &
            bind(c,name = 'mgis_bv_behaviour_get_internal_state_variable_type') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         implicit none
         integer t
         type(c_ptr), intent(in), value :: b
         integer(kind=c_size_t), intent(in), value :: n
         type(mgis_status) :: r
       end function behaviour_get_internal_state_variable_type_wrapper
    end interface
    integer, intent (out) :: t
    type(behaviour), intent(in) :: b
    integer, intent (in) :: n
    type(mgis_status) :: s
    integer(kind = c_size_t) :: nc
        if(.not. convert_to_c_index(nc, n)) then
       s = report_failure("invalid index")
       return
    end if
    s = behaviour_get_internal_state_variable_type_wrapper(t, b%ptr, nc)
  end function behaviour_get_internal_state_variable_type
  !
  function behaviour_get_internal_state_variable_offset(o, b, n) result(s)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function behaviour_get_internal_state_variable_offset_wrapper(o, b, n) &
            bind(c,name = 'mgis_bv_behaviour_get_internal_state_variable_offset') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t, c_char
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_size_t), intent(out) :: o
         type(c_ptr), intent(in), value :: b
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: r
       end function behaviour_get_internal_state_variable_offset_wrapper
    end interface
    integer, intent (out) :: o
    type(behaviour), intent(in) :: b
    character(len=*), intent(in) :: n
    type(mgis_status) :: s
    integer(kind=c_size_t) offset
    s = behaviour_get_internal_state_variable_offset_wrapper( &
         offset, b%ptr, convert_fortran_string(n))
    o = offset+1
  end function behaviour_get_internal_state_variable_offset
  !
  function behaviour_get_number_of_external_state_variables(n,b) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    implicit none
    interface
       function behaviour_get_number_of_external_state_variables_wrapper(l,b) &
            bind(c,name = 'mgis_bv_behaviour_get_number_of_external_state_variables') result(r)
         use, intrinsic :: iso_c_binding, only: c_size_t, c_ptr
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_size_t), intent(out) :: l
         type(c_ptr), intent(in), value :: b
         type(mgis_status) :: r
       end function behaviour_get_number_of_external_state_variables_wrapper
    end interface
    integer :: n
    type(behaviour), intent(in) :: b
    type(mgis_status) :: s
    integer(kind=c_size_t) :: ns
    s = behaviour_get_number_of_external_state_variables_wrapper(ns, b%ptr)
    n = ns
  end function behaviour_get_number_of_external_state_variables
  !
  function behaviour_get_external_state_variable_name(l, b, n) result(s)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_get_external_state_variable_name_wrapper(l, b, n) &
            bind(c,name = 'mgis_bv_behaviour_get_external_state_variable_name') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: l
         type(c_ptr), intent(in), value :: b
         integer(kind=c_size_t), intent(in), value :: n
         type(mgis_status) :: r
       end function behaviour_get_external_state_variable_name_wrapper
    end interface
    character(len=:), allocatable, intent(out) :: l
    type(behaviour), intent(in) :: b
    integer, intent (in) :: n
    type(mgis_status) :: s
    type(c_ptr) :: lp
    integer(kind = c_size_t) :: nc
    if(.not. convert_to_c_index(nc, n)) then
       s = report_failure("invalid index")
       return
    end if
    s = behaviour_get_external_state_variable_name_wrapper(lp, b%ptr, nc)
    if (s % exit_status == MGIS_SUCCESS) then
       l = convert_c_string(lp)
    else
       l = get_empty_string();
    end if
  end function behaviour_get_external_state_variable_name
  !
  function behaviour_get_external_state_variable_type(t, b, n) result(s)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function behaviour_get_external_state_variable_type_wrapper(t, b, n) &
            bind(c,name = 'mgis_bv_behaviour_get_external_state_variable_type') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         implicit none
         integer t
         type(c_ptr), intent(in), value :: b
         integer(kind=c_size_t), intent(in), value :: n
         type(mgis_status) :: r
       end function behaviour_get_external_state_variable_type_wrapper
    end interface
    integer, intent (out) :: t
    type(behaviour), intent(in) :: b
    integer, intent (in) :: n
    type(mgis_status) :: s
    integer(kind = c_size_t) :: nc
    if(.not. convert_to_c_index(nc, n)) then
       s = report_failure("invalid index")
       return
    end if
    s = behaviour_get_external_state_variable_type_wrapper(t, b%ptr, nc)
  end function behaviour_get_external_state_variable_type
  !
  function behaviour_has_bounds(r, b, n) result(s)
    use, intrinsic :: iso_c_binding, only: c_int
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_has_bounds_wrapper(r, b, n) &
            bind(c,name='mgis_bv_behaviour_has_bounds') result(s)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_int
         use mgis, only: mgis_status
         integer(kind=c_int), intent(out) :: r
         type(c_ptr), intent(in), value :: b
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: s
       end function behaviour_has_bounds_wrapper
    end interface
    logical, intent(out) :: r
    type(Behaviour), intent(in) :: b
    character(len=*), intent(in) :: n
    type(mgis_status) :: s
    integer(kind=c_int) :: rc
    s = behaviour_has_bounds_wrapper( &
         rc, b%ptr , convert_fortran_string(n))
    if( s % exit_status .eq. MGIS_SUCCESS) then
       r = rc .eq. 1
    else
       r = .false.
    end if
  end function behaviour_has_bounds
  !
  function behaviour_has_lower_bound(r, b, n) result(s)
    use, intrinsic :: iso_c_binding, only: c_int
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_has_lower_bound_wrapper(r, b, n) &
            bind(c,name='mgis_bv_behaviour_has_lower_bound') result(s)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_int
         use mgis, only: mgis_status
         integer(kind=c_int), intent(out) :: r
         type(c_ptr), intent(in), value :: b
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: s
       end function behaviour_has_lower_bound_wrapper
    end interface
    logical, intent(out) :: r
    type(Behaviour), intent(in) :: b
    character(len=*), intent(in) :: n
    type(mgis_status) :: s
    integer(kind=c_int) :: rc
    s = behaviour_has_lower_bound_wrapper( &
         rc, b%ptr , convert_fortran_string(n))
    if( s % exit_status .eq. MGIS_SUCCESS) then
       r = rc .eq. 1
    else
       r = .false.
    end if
  end function behaviour_has_lower_bound
  !
  function behaviour_get_lower_bound(v, b, n) result(s)
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_c_binding, only: c_long_double
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_get_lower_bound_wrapper(v, b, n) &
            bind(c,name='mgis_bv_behaviour_get_lower_bound') result(s)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_long_double
         use mgis, only: mgis_status
         real(kind=c_long_double), intent(out) :: v
         type(c_ptr), intent(in), value :: b
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: s
       end function behaviour_get_lower_bound_wrapper
    end interface
    real(kind=8), intent(out) :: v
    type(Behaviour), intent(in) :: b
    character(len=*), intent(in) :: n
    type(mgis_status) :: s
    real(kind=c_long_double) :: vc
    s = behaviour_get_lower_bound_wrapper( &
         vc, b%ptr , convert_fortran_string(n))
    if( s % exit_status .eq. MGIS_SUCCESS) then
       v = vc
    else
       v = ieee_value(v, ieee_quiet_nan)
    end if
  end function behaviour_get_lower_bound
  !
  function behaviour_has_upper_bound(r, b, n) result(s)
    use, intrinsic :: iso_c_binding, only: c_int
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_has_upper_bound_wrapper(r, b, n) &
            bind(c,name='mgis_bv_behaviour_has_upper_bound') result(s)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_int
         use mgis, only: mgis_status
         integer(kind=c_int), intent(out) :: r
         type(c_ptr), intent(in), value :: b
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: s
       end function behaviour_has_upper_bound_wrapper
    end interface
    logical, intent(out) :: r
    type(Behaviour), intent(in) :: b
    character(len=*), intent(in) :: n
    type(mgis_status) :: s
    integer(kind=c_int) :: rc
    s = behaviour_has_upper_bound_wrapper( &
         rc, b%ptr , convert_fortran_string(n))
    if( s % exit_status .eq. MGIS_SUCCESS) then
       r = rc .eq. 1
    else
       r = .false.
    end if
  end function behaviour_has_upper_bound
  !
  function behaviour_get_upper_bound(v, b, n) result(s)
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_c_binding, only: c_long_double
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_get_upper_bound_wrapper(v, b, n) &
            bind(c,name='mgis_bv_behaviour_get_upper_bound') result(s)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_long_double
         use mgis, only: mgis_status
         real(kind=c_long_double), intent(out) :: v
         type(c_ptr), intent(in), value :: b
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: s
       end function behaviour_get_upper_bound_wrapper
    end interface
    real(kind=8), intent(out) :: v
    type(Behaviour), intent(in) :: b
    character(len=*), intent(in) :: n
    type(mgis_status) :: s
    real(kind=c_long_double) :: vc
    s = behaviour_get_upper_bound_wrapper( &
         vc, b%ptr , convert_fortran_string(n))
    if( s % exit_status .eq. MGIS_SUCCESS) then
       v = vc
    else
       v = ieee_value(v, ieee_quiet_nan)
    end if
  end function behaviour_get_upper_bound
  !
  function behaviour_has_physical_bounds(r, b, n) result(s)
    use, intrinsic :: iso_c_binding, only: c_int
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_has_physical_bounds_wrapper(r, b, n) &
            bind(c,name='mgis_bv_behaviour_has_physical_bounds') result(s)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_int
         use mgis, only: mgis_status
         integer(kind=c_int), intent(out) :: r
         type(c_ptr), intent(in), value :: b
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: s
       end function behaviour_has_physical_bounds_wrapper
    end interface
    logical, intent(out) :: r
    type(Behaviour), intent(in) :: b
    character(len=*), intent(in) :: n
    type(mgis_status) :: s
    integer(kind=c_int) :: rc
    s = behaviour_has_physical_bounds_wrapper( &
         rc, b%ptr , convert_fortran_string(n))
    if( s % exit_status .eq. MGIS_SUCCESS) then
       r = rc .eq. 1
    else
       r = .false.
    end if
  end function behaviour_has_physical_bounds
  !
  function behaviour_has_lower_physical_bound(r, b, n) result(s)
    use, intrinsic :: iso_c_binding, only: c_int
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_has_lower_physical_bound_wrapper(r, b, n) &
            bind(c,name='mgis_bv_behaviour_has_lower_physical_bound') result(s)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_int
         use mgis, only: mgis_status
         integer(kind=c_int), intent(out) :: r
         type(c_ptr), intent(in), value :: b
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: s
       end function behaviour_has_lower_physical_bound_wrapper
    end interface
    logical, intent(out) :: r
    type(Behaviour), intent(in) :: b
    character(len=*), intent(in) :: n
    type(mgis_status) :: s
    integer(kind=c_int) :: rc
    s = behaviour_has_lower_physical_bound_wrapper( &
         rc, b%ptr , convert_fortran_string(n))
    if( s % exit_status .eq. MGIS_SUCCESS) then
       r = rc .eq. 1
    else
       r = .false.
    end if
  end function behaviour_has_lower_physical_bound
  !
  function behaviour_get_lower_physical_bound(v, b, n) result(s)
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_c_binding, only: c_long_double
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_get_lower_physical_bound_wrapper(v, b, n) &
            bind(c,name='mgis_bv_behaviour_get_lower_physical_bound') result(s)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_long_double
         use mgis, only: mgis_status
         real(kind=c_long_double), intent(out) :: v
         type(c_ptr), intent(in), value :: b
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: s
       end function behaviour_get_lower_physical_bound_wrapper
    end interface
    real(kind=8), intent(out) :: v
    type(Behaviour), intent(in) :: b
    character(len=*), intent(in) :: n
    type(mgis_status) :: s
    real(kind=c_long_double) :: vc
    s = behaviour_get_lower_physical_bound_wrapper( &
         vc, b%ptr , convert_fortran_string(n))
    if( s % exit_status .eq. MGIS_SUCCESS) then
       v = vc
    else
       v = ieee_value(v, ieee_quiet_nan)
    end if
  end function behaviour_get_lower_physical_bound
  !
  function behaviour_has_upper_physical_bound(r, b, n) result(s)
    use, intrinsic :: iso_c_binding, only: c_int
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    interface
       function behaviour_has_upper_physical_bound_wrapper(r, b, n) &
            bind(c,name='mgis_bv_behaviour_has_upper_physical_bound') result(s)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_int
         use mgis, only: mgis_status
         integer(kind=c_int), intent(out) :: r
         type(c_ptr), intent(in), value :: b
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: s
       end function behaviour_has_upper_physical_bound_wrapper
    end interface
    logical, intent(out) :: r
    type(Behaviour), intent(in) :: b
    character(len=*), intent(in) :: n
    type(mgis_status) :: s
    integer(kind=c_int) :: rc
    s = behaviour_has_upper_physical_bound_wrapper( &
         rc, b%ptr , convert_fortran_string(n))
    if( s % exit_status .eq. MGIS_SUCCESS) then
       r = rc .eq. 1
    else
       r = .false.
    end if
  end function behaviour_has_upper_physical_bound
  !
  function behaviour_get_upper_physical_bound(v, b, n) result(s)
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_c_binding, only: c_long_double
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_get_upper_physical_bound_wrapper(v, b, n) &
            bind(c,name='mgis_bv_behaviour_get_upper_physical_bound') result(s)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_char, c_long_double
         use mgis, only: mgis_status
         real(kind=c_long_double), intent(out) :: v
         type(c_ptr), intent(in), value :: b
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: s
       end function behaviour_get_upper_physical_bound_wrapper
    end interface
    real(kind=8), intent(out) :: v
    type(Behaviour), intent(in) :: b
    character(len=*), intent(in) :: n
    type(mgis_status) :: s
    real(kind=c_long_double) :: vc
    s = behaviour_get_upper_physical_bound_wrapper( &
         vc, b%ptr , convert_fortran_string(n))
    if( s % exit_status .eq. MGIS_SUCCESS) then
       v = vc
    else
       v = ieee_value(v, ieee_quiet_nan)
    end if
  end function behaviour_get_upper_physical_bound
  ! free behaviour
  function free_behaviour(ptr) result(r)
    use, intrinsic :: iso_c_binding, only: c_associated
    use mgis
    implicit none
    interface
       function free_behaviour_wrapper(ptr) bind(c, name='mgis_bv_free_behaviour') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis
         implicit none
         type(c_ptr), intent(inout) :: ptr
         type(mgis_status) :: r
       end function free_behaviour_wrapper
    end interface
    type(Behaviour), intent(inout) :: ptr
    type(mgis_status) :: r
    if (c_associated(ptr%ptr)) then
       r = free_behaviour_wrapper(ptr%ptr)
    end if
  end function free_behaviour
  !
  function allocate_behaviour_data(d,b) result(s)
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    implicit none
    interface
       function allocate_behaviour_data_wrapper(d,b) bind(c,name = 'mgis_bv_allocate_behaviour_data') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: d
         type(c_ptr), intent(in), value :: b
         type(mgis_status) :: r
       end function allocate_behaviour_data_wrapper
    end interface
    type(BehaviourData), intent(out) :: d
    type(Behaviour), intent(in) :: b
    type(mgis_status) :: s
    s = allocate_behaviour_data_wrapper(d%ptr, b%ptr)
  end function allocate_behaviour_data
  !
  function behaviour_data_get_behaviour(b, d) result(s)
    use mgis, only: mgis_status
    implicit none
    interface
       function behaviour_data_get_behaviour_wrapper(b, d) &
            bind(c,name = 'mgis_bv_behaviour_data_get_behaviour') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: b
         type(c_ptr), intent(in), value :: d
         type(mgis_status) :: r
       end function behaviour_data_get_behaviour_wrapper
    end interface
    type(Behaviour), intent(out) :: b
    type(BehaviourData), intent(in) :: d
    type(mgis_status) :: s
    s = behaviour_data_get_behaviour_wrapper(b%ptr, d%ptr)
  end function behaviour_data_get_behaviour
  !
  function behaviour_data_get_state_0(s0,d) result(s)
    use mgis, only: mgis_status
    implicit none
    interface
       function behaviour_data_get_state_0_wrapper(s0, d) &
            bind(c,name = 'mgis_bv_behaviour_data_get_state_0') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: s0
         type(c_ptr), intent(in), value :: d
         type(mgis_status) :: r
       end function behaviour_data_get_state_0_wrapper
    end interface
    type(State), intent(out) :: s0
    type(BehaviourData), intent(in) :: d
    type(mgis_status) :: s
    s = behaviour_data_get_state_0_wrapper(s0%ptr, d%ptr)
  end function behaviour_data_get_state_0
  !
  function behaviour_data_get_state_1(s1,d) result(s)
    use mgis, only: mgis_status
    implicit none
    interface
       function behaviour_data_get_state_1_wrapper(s1, d) &
            bind(c,name = 'mgis_bv_behaviour_data_get_state_1') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: s1
         type(c_ptr), intent(in), value :: d
         type(mgis_status) :: r
       end function behaviour_data_get_state_1_wrapper
    end interface
    type(State), intent(out) :: s1
    type(BehaviourData), intent(in) :: d
    type(mgis_status) :: s
    s = behaviour_data_get_state_1_wrapper(s1%ptr, d%ptr)
  end function behaviour_data_get_state_1
  !
  function behaviour_data_set_time_increment(d, dt) result(s)
    use mgis, only: mgis_status
    implicit none
    interface
       function behaviour_data_set_time_increment_wrapper(d, v) &
            bind(c,name = 'mgis_bv_behaviour_data_set_time_increment') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_double
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(in), value :: d
         real(kind=c_double), intent(in), value :: v
         type(mgis_status) :: r
       end function behaviour_data_set_time_increment_wrapper
    end interface
    type(BehaviourData), intent(in) :: d
    real(kind = 8) :: dt
    type(mgis_status) :: s
    s = behaviour_data_set_time_increment_wrapper(d%ptr, dt)
  end function behaviour_data_set_time_increment
  !
  function behaviour_data_get_time_step_scaling_factor(rdt, d) result(s)
    use mgis, only: mgis_status
    implicit none
    interface
       function behaviour_data_get_time_step_scaling_factor_wrapper(v, d) &
            bind(c,name = 'mgis_bv_behaviour_data_get_time_step_scaling_factor')&
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_double
         use mgis, only: mgis_status
         implicit none
         real(kind=c_double), intent(out) :: v
         type(c_ptr), intent(in), value :: d
         type(mgis_status) :: r
       end function behaviour_data_get_time_step_scaling_factor_wrapper
    end interface
    real(kind = 8) :: rdt
    type(BehaviourData), intent(in) :: d
    type(mgis_status) :: s
    s = behaviour_data_get_time_step_scaling_factor_wrapper(rdt, d%ptr)
  end function behaviour_data_get_time_step_scaling_factor
  ! \brief return the transpose of the tangent operator
  function behaviour_data_get_tangent_operator(K, d) result(s)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function behaviour_data_get_tangent_operator_wrapper(K, d) &
            bind(c,name = 'mgis_bv_behaviour_data_get_tangent_operator')&
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: K
         type(c_ptr), intent(in), value :: d
         type(mgis_status) :: r
       end function behaviour_data_get_tangent_operator_wrapper
    end interface
    real(kind=8), dimension(:,:), pointer, intent(out) :: K
    type(BehaviourData), intent(in) :: d
    type(mgis_status) :: s
    type(c_ptr) :: p
    type(Behaviour) :: b
    integer gs
    integer ths
    nullify(K)
    s = behaviour_data_get_behaviour(b, d)
    if( s%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    s = behaviour_get_gradients_size(gs, b)
    if( s%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    s = behaviour_get_thermodynamic_forces_size(ths, b)
    if( s%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    s = behaviour_data_get_tangent_operator_wrapper(p, d%ptr)
    if( s%exit_status .eq. MGIS_SUCCESS) then
       call c_f_pointer(p,K,[ths,gs])
    end if
  end function behaviour_data_get_tangent_operator
  !
  function update_behaviour_data(d) result(r)
    use mgis
    implicit none
    interface
       function update_behaviour_data_wrapper(ptr) &
            bind(c, name='mgis_bv_update_behaviour_data') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis
         implicit none
         type(c_ptr), intent(in), value :: ptr
         type(mgis_status) :: r
       end function update_behaviour_data_wrapper
    end interface
    type(BehaviourData), intent(in) :: d
    type(mgis_status) :: r
    r = update_behaviour_data_wrapper(d%ptr)
  end function update_behaviour_data
  !
  function revert_behaviour_data(d) result(r)
    use mgis
    implicit none
    interface
       function revert_behaviour_data_wrapper(ptr) &
            bind(c, name='mgis_bv_revert_behaviour_data') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis
         implicit none
         type(c_ptr), intent(in), value :: ptr
         type(mgis_status) :: r
       end function revert_behaviour_data_wrapper
    end interface
    type(BehaviourData), intent(in) :: d
    type(mgis_status) :: r
    r = revert_behaviour_data_wrapper(d%ptr)
  end function revert_behaviour_data
  !
  function free_behaviour_data(ptr) result(r)
    use, intrinsic :: iso_c_binding, only: c_associated
    use mgis
    implicit none
    interface
       function free_behaviour_data_wrapper(ptr) bind(c, name='mgis_bv_free_behaviour_data') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis
         implicit none
         type(c_ptr), intent(inout) :: ptr
         type(mgis_status) :: r
       end function free_behaviour_data_wrapper
    end interface
    type(BehaviourData), intent(inout) :: ptr
    type(mgis_status) :: r
    if (c_associated(ptr%ptr)) then
       r = free_behaviour_data_wrapper(ptr%ptr)
    end if
  end function free_behaviour_data
  !
  function state_set_gradient_by_name(s, n, v) result(r)
    use, intrinsic :: iso_c_binding, only: c_loc
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    implicit none
    interface
       function state_set_gradient_by_name_wrapper(s, n, v) &
            bind(c,name = 'mgis_bv_state_set_gradient_by_name') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_char
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(in),value :: s
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(c_ptr), intent(in),value :: v
         type(mgis_status) :: r
       end function state_set_gradient_by_name_wrapper
    end interface
    type(State), intent(in) :: s
    character(len=*), intent(in) :: n
    real(kind=8), dimension(*), target :: v
    type(mgis_status) :: r
    r = state_set_gradient_by_name_wrapper( &
         s%ptr, convert_fortran_string(n), c_loc(v(1)))
  end function state_set_gradient_by_name
  !
  function state_set_material_property_by_name(s, n, v) result(r)
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    implicit none
    interface
       function state_set_material_property_by_name_wrapper(s1, n, v) &
            bind(c,name = 'mgis_bv_state_set_material_property_by_name') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_char
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(in),value :: s1
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         real(kind=c_double), intent(in), value :: v
         type(mgis_status) :: r
       end function state_set_material_property_by_name_wrapper
    end interface
    type(State), intent(in) :: s
    character(len=*), intent(in) :: n
    real(kind=8), intent(in) :: v
    type(mgis_status) :: r
    r = state_set_material_property_by_name_wrapper(&
         s%ptr, convert_fortran_string(n), v)
  end function state_set_material_property_by_name
  ! \note the C function returns a pointer to the variable, the
  !       fortran funtion returns the value
  function state_get_material_property_by_name(v, s, n) result(r)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function state_get_material_property_by_name_wrapper(v, s1, n) &
            bind(c,name = 'mgis_bv_state_get_material_property_by_name') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_char
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: v
         type(c_ptr), intent(in),value :: s1
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: r
       end function state_get_material_property_by_name_wrapper
    end interface
    real(kind=8), intent(out) :: v
    type(State), intent(in) :: s
    character(len=*), intent(in) :: n
    type(mgis_status) :: r
    real(kind=8), pointer :: fptr => null()
    type(c_ptr) :: p
    r = state_get_material_property_by_name_wrapper(&
         p, s%ptr, convert_fortran_string(n))
    if (r % exit_status == MGIS_SUCCESS) then
       call c_f_pointer(p,fptr)
       v = fptr
    else
       v = ieee_value(v, ieee_quiet_nan)
    endif
  end function state_get_material_property_by_name
  ! \note the C function returns a pointer to the variable, the
  !       fortran funtion returns the value
  function state_get_internal_state_variable_by_offset(v, s, o) result(r)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_size_t
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure, MGIS_SUCCESS
    implicit none
    interface
       function state_get_internal_state_variable_by_offset_wrapper(v, s1, o) &
            bind(c,name = 'mgis_bv_state_get_internal_state_variable_by_offset') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: v
         type(c_ptr), intent(in),value :: s1
         integer(kind=c_size_t), intent(in), value :: o
         type(mgis_status) :: r
       end function state_get_internal_state_variable_by_offset_wrapper
    end interface
    real(kind=8), intent(out) :: v
    type(State), intent(in) :: s
    integer, intent(in) :: o
    type(mgis_status) :: r
    real(kind=8), pointer :: fptr => null()
    type(c_ptr) :: p
    integer(kind=c_size_t) offset
    if(.not. convert_to_c_index(offset, o)) then
       r = report_failure("invalid index")
       return
    end if
    r = state_get_internal_state_variable_by_offset_wrapper(&
         p, s%ptr, offset)
    if (r % exit_status == MGIS_SUCCESS) then
       call c_f_pointer(p,fptr)
       v = fptr
    else
       v = ieee_value(v, ieee_quiet_nan)
    endif
  end function state_get_internal_state_variable_by_offset
  !
  function state_set_external_state_variable_by_name(s, n, v) result(r)
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    implicit none
    interface
       function state_set_external_state_variable_by_name_wrapper(s1, n, v) &
            bind(c,name = 'mgis_bv_state_set_external_state_variable_by_name') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_char
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(in),value :: s1
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         real(kind=c_double), intent(in), value :: v
         type(mgis_status) :: r
       end function state_set_external_state_variable_by_name_wrapper
    end interface
    type(State), intent(in) :: s
    character(len=*), intent(in) :: n
    real(kind=8), intent(in) :: v
    type(mgis_status) :: r
    r = state_set_external_state_variable_by_name_wrapper(&
         s%ptr, convert_fortran_string(n), v)
  end function state_set_external_state_variable_by_name
  ! \note the C function returns a pointer to the variable, the
  !       fortran funtion returns the value
  function state_get_external_state_variable_by_name(v, s, n) result(r)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function state_get_external_state_variable_by_name_wrapper(v, s1, n) &
            bind(c,name = 'mgis_bv_state_get_external_state_variable_by_name') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_char
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: v
         type(c_ptr), intent(in),value :: s1
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: r
       end function state_get_external_state_variable_by_name_wrapper
    end interface
    real(kind=8), intent(out) :: v
    type(State), intent(in) :: s
    character(len=*), intent(in) :: n
    type(mgis_status) :: r
    real(kind=8), pointer :: fptr => null()
    type(c_ptr) :: p
    r = state_get_external_state_variable_by_name_wrapper(&
         p, s%ptr, convert_fortran_string(n))
    if (r % exit_status == MGIS_SUCCESS) then
       call c_f_pointer(p,fptr)
       v = fptr
    else
       v = ieee_value(v, ieee_quiet_nan)
    endif
  end function state_get_external_state_variable_by_name
  !
  function msm_initializer_bind_gradients(st,a,n) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t, c_loc
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function msm_initializer_bind_gradients_wrapper(stc,ac,nc) &
            bind(c,name = 'mgis_bv_material_state_manager_initializer_bind_gradients') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(in), value :: stc
         type(c_ptr), intent(in), value :: ac
         integer(kind=c_size_t), intent(in), value :: nc
         type(mgis_status) :: r
       end function msm_initializer_bind_gradients_wrapper
    end interface
    type(MaterialStateManagerInitializer), intent(in) :: st
    real(kind=8), dimension(*), target, intent(out) :: a
    integer, intent(in) :: n
    type(mgis_status) :: s
    type(c_ptr) a_ptr
    integer(kind=c_size_t) nc
    nc=n
    a_ptr = c_loc(a)
    s = msm_initializer_bind_gradients_wrapper(st%ptr,a_ptr,nc)
  end function msm_initializer_bind_gradients
  !
  function msm_initializer_bind_thermodynamic_forces(st,a,n) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t, c_loc
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function msm_initializer_bind_thermodynamic_forces_wrapper(stc,ac,nc) &
            bind(c,name = 'mgis_bv_material_state_manager_initializer_bind_thermodynamic_forces') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(in), value :: stc
         type(c_ptr), intent(in), value :: ac
         integer(kind=c_size_t), intent(in), value :: nc
         type(mgis_status) :: r
       end function msm_initializer_bind_thermodynamic_forces_wrapper
    end interface
    type(MaterialStateManagerInitializer), intent(in) :: st
    real(kind=8), dimension(*), target, intent(out) :: a
    integer, intent(in) :: n
    type(mgis_status) :: s
    type(c_ptr) a_ptr
    integer(kind=c_size_t) nc
    nc=n
    a_ptr = c_loc(a)
    s = msm_initializer_bind_thermodynamic_forces_wrapper(st%ptr,a_ptr,nc)
  end function msm_initializer_bind_thermodynamic_forces
  !
  function msm_initializer_bind_internal_state_variables(st,a,n) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t, c_loc
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function msm_initializer_bind_internal_state_variables_wrapper(stc,ac,nc) &
            bind(c,name = 'mgis_bv_material_state_manager_initializer_bind_internal_state_variables') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(in), value :: stc
         type(c_ptr), intent(in), value :: ac
         integer(kind=c_size_t), intent(in), value :: nc
         type(mgis_status) :: r
       end function msm_initializer_bind_internal_state_variables_wrapper
    end interface
    type(MaterialStateManagerInitializer), intent(in) :: st
    real(kind=8), dimension(*), target, intent(out) :: a
    integer, intent(in) :: n
    type(mgis_status) :: s
    type(c_ptr) a_ptr
    integer(kind=c_size_t) nc
    nc=n
    a_ptr = c_loc(a)
    s = msm_initializer_bind_internal_state_variables_wrapper(st%ptr,a_ptr,nc)
  end function msm_initializer_bind_internal_state_variables
  !
  function msm_initializer_bind_stored_energies(st,a,n) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t, c_loc
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function msm_initializer_bind_stored_energies_wrapper(stc,ac,nc) &
            bind(c,name = 'mgis_bv_material_state_manager_initializer_bind_stored_energies') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(in), value :: stc
         type(c_ptr), intent(in), value :: ac
         integer(kind=c_size_t), intent(in), value :: nc
         type(mgis_status) :: r
       end function msm_initializer_bind_stored_energies_wrapper
    end interface
    type(MaterialStateManagerInitializer), intent(in) :: st
    real(kind=8), dimension(*), target, intent(out) :: a
    integer, intent(in) :: n
    type(mgis_status) :: s
    type(c_ptr) a_ptr
    integer(kind=c_size_t) nc
    nc=n
    a_ptr = c_loc(a)
    s = msm_initializer_bind_stored_energies_wrapper(st%ptr,a_ptr,nc)
  end function msm_initializer_bind_stored_energies
  !
  function msm_initializer_bind_dissipated_energies(st,a,n) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t, c_loc
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function msm_initializer_bind_dissipated_energies_wrapper(stc,ac,nc) &
            bind(c,name = 'mgis_bv_material_state_manager_initializer_bind_dissipated_energies') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(in), value :: stc
         type(c_ptr), intent(in), value :: ac
         integer(kind=c_size_t), intent(in), value :: nc
         type(mgis_status) :: r
       end function msm_initializer_bind_dissipated_energies_wrapper
    end interface
    type(MaterialStateManagerInitializer), intent(in) :: st
    real(kind=8), dimension(*), target, intent(out) :: a
    integer, intent(in) :: n
    type(mgis_status) :: s
    type(c_ptr) a_ptr
    integer(kind=c_size_t) nc
    nc=n
    a_ptr = c_loc(a)
    s = msm_initializer_bind_dissipated_energies_wrapper(st%ptr,a_ptr,nc)
  end function msm_initializer_bind_dissipated_energies
  !
  function create_material_data_manager_initializer(d) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function create_material_data_manager_initializer_wrapper(d) &
            bind(c,name = 'mgis_bv_create_material_data_manager_initializer') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: d
         type(mgis_status) :: r
       end function create_material_data_manager_initializer_wrapper
    end interface
    type(MaterialDataManagerInitializer), intent(out) :: d
    type(mgis_status) :: s
    s = create_material_data_manager_initializer_wrapper(d%ptr)
  end function create_material_data_manager_initializer
  !
  function material_data_manager_initializer_get_state_0_initializer(s0,d) result(s)
    use mgis, only: mgis_status
    implicit none
    interface
       function mdm_initializer_get_state_0_initializer_wrapper(s0c,dc) &
            bind(c,name = 'mgis_bv_material_data_manager_initializer_get_state_0_initializer') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: s0c
         type(c_ptr), intent(in), value :: dc
         type(mgis_status) :: r
       end function mdm_initializer_get_state_0_initializer_wrapper
    end interface
    type(MaterialDataManagerInitializer), intent(in) :: d
    type(MaterialStateManagerInitializer), intent(out) :: s0
    type(mgis_status) :: s
    s = mdm_initializer_get_state_0_initializer_wrapper(s0%ptr, d%ptr)
  end function material_data_manager_initializer_get_state_0_initializer
  !
  function material_data_manager_initializer_get_state_1_initializer(s1,d) result(s)
    use mgis, only: mgis_status
    implicit none
    interface
       function mdm_initializer_get_state_1_initializer_wrapper(s1c,dc) &
            bind(c,name = 'mgis_bv_material_data_manager_initializer_get_state_1_initializer') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: s1c
         type(c_ptr), intent(in), value :: dc
         type(mgis_status) :: r
       end function mdm_initializer_get_state_1_initializer_wrapper
    end interface
    type(MaterialDataManagerInitializer), intent(in) :: d
    type(MaterialStateManagerInitializer), intent(out) :: s1
    type(mgis_status) :: s
    s = mdm_initializer_get_state_1_initializer_wrapper(s1%ptr, d%ptr)
  end function material_data_manager_initializer_get_state_1_initializer 
  !
  function material_data_manager_initializer_bind_tangent_operator(d,K,n) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t, c_loc
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function material_data_manager_initializer_bind_tangent_operator_wrapper(d,K,s) &
            bind(c,name = 'mgis_bv_material_data_manager_initializer_bind_tangent_operator') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(in), value :: d
         type(c_ptr), intent(in), value :: K
         integer(kind=c_size_t), intent(in), value :: s
         type(mgis_status) :: r
       end function material_data_manager_initializer_bind_tangent_operator_wrapper
    end interface
    type(MaterialDataManagerInitializer), intent(in) :: d
    real(kind=8), dimension(:,:), target, intent(out) :: K
    integer, intent(in) :: n
    type(mgis_status) :: s
    type(c_ptr) K_ptr
    integer(kind=c_size_t) nc
    nc=n
    K_ptr = c_loc(K)
    s = material_data_manager_initializer_bind_tangent_operator_wrapper(d%ptr,K_ptr,nc)
  end function material_data_manager_initializer_bind_tangent_operator
  !
  function free_material_data_manager_initializer(ptr) result(r)
    use, intrinsic :: iso_c_binding, only: c_associated
    use mgis
    implicit none
    interface
       function free_material_data_manager_initializer_wrapper(ptr) &
            bind(c, name='mgis_bv_free_material_data_manager_initializer') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis
         implicit none
         type(c_ptr), intent(inout) :: ptr
         type(mgis_status) :: r
       end function free_material_data_manager_initializer_wrapper
    end interface
    type(MaterialDataManagerInitializer), intent(inout) :: ptr
    type(mgis_status) :: r
    if (c_associated(ptr%ptr)) then
       r = free_material_data_manager_initializer_wrapper(ptr%ptr)
    end if
  end function free_material_data_manager_initializer
  !
  function create_material_data_manager(d, b, n) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function create_material_data_manager_wrapper(d, b, n) &
            bind(c,name = 'mgis_bv_create_material_data_manager') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: d
         type(c_ptr), intent(in), value :: b
         integer(kind=c_size_t), intent(in), value :: n
         type(mgis_status) :: r
       end function create_material_data_manager_wrapper
    end interface
    type(MaterialDataManager), intent(out) :: d
    type(Behaviour), intent(in) :: b
    integer, intent(in) :: n
    type(mgis_status) :: s
    integer(kind=c_size_t) nc
    if (n.lt.1) then
       s = report_failure('invalid number of integration points')
       return
    end if
    nc = n
    s = create_material_data_manager_wrapper(d%ptr, b%ptr, nc)
  end function create_material_data_manager
  !
  function create_material_data_manager_with_initializer(d, b, i, n) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t
    use mgis_fortran_utilities
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function create_material_data_manager_with_initializer_wrapper(d, b, n, i) &
            bind(c,name = 'mgis_bv_create_material_data_manager_with_initializer') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: d
         type(c_ptr), intent(in), value :: b
         type(c_ptr), intent(in), value :: i
         integer(kind=c_size_t), intent(in), value :: n
         type(mgis_status) :: r
       end function create_material_data_manager_with_initializer_wrapper
    end interface
    type(MaterialDataManager), intent(out) :: d
    type(Behaviour), intent(in) :: b
    type(MaterialDataManagerInitializer), intent(in) :: i
    integer, intent(in) :: n
    type(mgis_status) :: s
    integer(kind=c_size_t) nc
    if (n.lt.1) then
       s = report_failure('invalid number of integration points')
       return
    end if
    nc = n
    s = create_material_data_manager_with_initializer_wrapper(d%ptr, b%ptr, nc, i%ptr)
  end function create_material_data_manager_with_initializer
  !
  function material_data_manager_get_state_0(s0,d) result(s)
    use mgis, only: mgis_status
    implicit none
    interface
       function material_data_manager_get_state_0_wrapper(s0, d) &
            bind(c,name = 'mgis_bv_material_data_manager_get_state_0') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: s0
         type(c_ptr), intent(in), value :: d
         type(mgis_status) :: r
       end function material_data_manager_get_state_0_wrapper
    end interface
    type(MaterialStateManager), intent(out) :: s0
    type(MaterialDataManager), intent(in) :: d
    type(mgis_status) :: s
    s = material_data_manager_get_state_0_wrapper(s0%ptr, d%ptr)
  end function material_data_manager_get_state_0
  !
  function material_data_manager_get_state_1(s1,d) result(s)
    use mgis, only: mgis_status
    implicit none
    interface
       function material_data_manager_get_state_1_wrapper(s1, d) &
            bind(c,name = 'mgis_bv_material_data_manager_get_state_1') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: s1
         type(c_ptr), intent(in), value :: d
         type(mgis_status) :: r
       end function material_data_manager_get_state_1_wrapper
    end interface
    type(MaterialStateManager), intent(out) :: s1
    type(MaterialDataManager), intent(in) :: d
    type(mgis_status) :: s
    s = material_data_manager_get_state_1_wrapper(s1%ptr, d%ptr)
  end function material_data_manager_get_state_1
  !
  function material_state_manager_get_number_of_integration_points(n, s) &
       result(r)
    use, intrinsic :: iso_c_binding, only: c_size_t
    use mgis, only: mgis_status, MGIS_SUCCESS
    use mgis_fortran_utilities
    implicit none
    interface
       function msm_get_number_of_integration_points_wrapper(nig, s) &
            bind(c,name = 'mgis_bv_material_state_manager_get_number_of_integration_points') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         integer(kind=c_size_t), intent(out) :: nig
         type(c_ptr), intent(in),value :: s
         type(mgis_status) :: r
       end function msm_get_number_of_integration_points_wrapper
    end interface
    integer, intent(out) :: n
    type(MaterialStateManager), intent(in) :: s
    type(mgis_status) :: r
    integer(kind=c_size_t) nig
    r = msm_get_number_of_integration_points_wrapper(nig, s %ptr)
    if( r%exit_status .eq. MGIS_SUCCESS) then
       n = nig
    else
       n = -1
    end if
  end function material_state_manager_get_number_of_integration_points
  !
  function material_state_manager_get_gradients_stride(g_stride, s) &
       result(r)
    use, intrinsic :: iso_c_binding, only: c_size_t
    use mgis, only: mgis_status, MGIS_SUCCESS
    use mgis_fortran_utilities
    implicit none
    interface
       function material_state_manager_get_gradients_stride_wrapper(gs, s) &
            bind(c,name = 'mgis_bv_material_state_manager_get_gradients_stride') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         integer(kind=c_size_t), intent(out) :: gs
         type(c_ptr), intent(in),value :: s
         type(mgis_status) :: r
       end function material_state_manager_get_gradients_stride_wrapper
    end interface
    integer, intent(out) :: g_stride
    type(MaterialStateManager), intent(in) :: s
    type(mgis_status) :: r
    integer(kind=c_size_t) gs
    r = material_state_manager_get_gradients_stride_wrapper(gs, s %ptr)
    if( r%exit_status .eq. MGIS_SUCCESS) then
       g_stride = gs
    else
       g_stride = -1
    end if
  end function material_state_manager_get_gradients_stride
  !
  function material_state_manager_get_gradients(g, s) &
       result(r)
    use, intrinsic :: iso_c_binding, only: c_size_t, c_f_pointer
    use mgis, only: mgis_status, MGIS_SUCCESS
    use mgis_fortran_utilities
    implicit none
    interface
       function material_state_manager_get_gradients_wrapper(g_ptr, s) &
            bind(c,name = 'mgis_bv_material_state_manager_get_gradients') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         type(c_ptr), intent(out) :: g_ptr
         type(c_ptr), intent(in),value :: s
         type(mgis_status) :: r
       end function material_state_manager_get_gradients_wrapper
    end interface
    real(kind=8), dimension(:,:), pointer, intent(out) :: g
    type(MaterialStateManager), intent(in) :: s
    type(mgis_status) :: r
    type(c_ptr) g_ptr
    integer gs
    integer n
    nullify(g)
    r = material_state_manager_get_gradients_wrapper(g_ptr, s %ptr)
    if( r%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    r = material_state_manager_get_gradients_stride(gs, s)
    if( r%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    r = material_state_manager_get_number_of_integration_points(n, s)
    if( r%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    call c_f_pointer(g_ptr, g, [gs, n])
  end function material_state_manager_get_gradients
  !
  function material_state_manager_get_thermodynamic_forces_stride(tf_stride, s) &
       result(r)
    use, intrinsic :: iso_c_binding, only: c_size_t
    use mgis, only: mgis_status, MGIS_SUCCESS
    use mgis_fortran_utilities
    implicit none
    interface
       function msm_get_thermodynamic_forces_stride_wrapper(tfs, s) &
            bind(c,name = 'mgis_bv_material_state_manager_get_thermodynamic_forces_stride') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         integer(kind=c_size_t), intent(out) :: tfs
         type(c_ptr), intent(in),value :: s
         type(mgis_status) :: r
       end function msm_get_thermodynamic_forces_stride_wrapper
    end interface
    integer, intent(out) :: tf_stride
    type(MaterialStateManager), intent(in) :: s
    type(mgis_status) :: r
    integer(kind=c_size_t) tfs
    r = msm_get_thermodynamic_forces_stride_wrapper(tfs, s %ptr)
    if( r%exit_status .eq. MGIS_SUCCESS) then
       tf_stride = tfs
    else
       tf_stride = -1
    end if
  end function material_state_manager_get_thermodynamic_forces_stride
  !
  function material_state_manager_get_thermodynamic_forces(tf, s) &
       result(r)
    use, intrinsic :: iso_c_binding, only: c_size_t, c_f_pointer
    use mgis, only: mgis_status, MGIS_SUCCESS
    use mgis_fortran_utilities
    implicit none
    interface
       function msm_get_thermodynamic_forces_wrapper(tf_ptr, s) &
            bind(c,name = 'mgis_bv_material_state_manager_get_thermodynamic_forces') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         type(c_ptr), intent(out) :: tf_ptr
         type(c_ptr), intent(in),value :: s
         type(mgis_status) :: r
       end function msm_get_thermodynamic_forces_wrapper
    end interface
    real(kind=8), dimension(:,:), pointer, intent(out) :: tf
    type(MaterialStateManager), intent(in) :: s
    type(mgis_status) :: r
    type(c_ptr) tf_ptr
    integer tfs
    integer n
    nullify(tf)
    r = msm_get_thermodynamic_forces_wrapper(tf_ptr, s %ptr)
    if( r%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    r = material_state_manager_get_thermodynamic_forces_stride(tfs, s)
    if( r%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    r = material_state_manager_get_number_of_integration_points(n, s)
    if( r%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    call c_f_pointer(tf_ptr, tf, [tfs, n])
  end function material_state_manager_get_thermodynamic_forces
  !
  function material_state_manager_set_uniform_material_property(s, n, v) &
       result(r)
    use mgis, only: mgis_status
    use mgis_fortran_utilities
    implicit none
    interface
       function msm_set_uniform_material_property_wrapper(s, n, v) &
            bind(c,name = 'mgis_bv_material_state_manager_set_uniform_material_property') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_char
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(in),value :: s
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         real(kind=c_double), intent(in), value :: v
         type(mgis_status) :: r
       end function msm_set_uniform_material_property_wrapper
    end interface
    type(MaterialStateManager), intent(in) :: s
    character(len=*), intent(in) :: n
    real(kind=8),     intent(in) :: v
    type(mgis_status) :: r
    r = msm_set_uniform_material_property_wrapper(s %ptr, &
         convert_fortran_string(n), v)
  end function material_state_manager_set_uniform_material_property
  !
  function material_state_manager_is_material_property_defined(b, s, n)&
       result(r)
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function msm_is_material_property_defined(b, s, n) &
            bind(c,name = 'mgis_bv_material_state_manager_is_material_property_defined') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_char
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_int), intent(out) :: b
         type(c_ptr), intent(in),value :: s
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: r
       end function msm_is_material_property_defined
    end interface
    logical, intent(out) :: b
    type(MaterialStateManager), intent(in) :: s
    character(len=*), intent(in) :: n
    type(mgis_status) :: r
    integer bc
    r = msm_is_material_property_defined(&
         bc, s%ptr, convert_fortran_string(n))
    if( r % exit_status .eq. MGIS_SUCCESS) then
       b = bc .eq. 1
    else
       b = .false.
    end if
  end function material_state_manager_is_material_property_defined
  !
  function material_state_manager_is_material_property_uniform(b, s, n)&
       result(r)
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function msm_is_material_property_uniform(b, s, n) &
            bind(c,name = 'mgis_bv_material_state_manager_is_material_property_uniform') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_char
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_int), intent(out) :: b
         type(c_ptr), intent(in),value :: s
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: r
       end function msm_is_material_property_uniform
    end interface
    logical, intent(out) :: b
    type(MaterialStateManager), intent(in) :: s
    character(len=*), intent(in) :: n
    type(mgis_status) :: r
    integer bc
    r = msm_is_material_property_uniform(&
         bc, s%ptr, convert_fortran_string(n))
    if( r % exit_status .eq. MGIS_SUCCESS) then
       b = bc .eq. 1
    else
       b = .false.
    end if
  end function material_state_manager_is_material_property_uniform
  !
  function material_state_manager_get_uniform_material_property(v, s, n)&
       result(r)
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function msm_get_uniform_material_property(v, s, n) &
            bind(c,name = 'mgis_bv_material_state_manager_get_uniform_material_property') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_char
         use mgis, only: mgis_status
         implicit none
         real(kind=c_double), intent(out) :: v
         type(c_ptr), intent(in),value :: s
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: r
       end function msm_get_uniform_material_property
    end interface
    real(kind=8), intent(out) :: v
    type(MaterialStateManager), intent(in) :: s
    character(len=*), intent(in) :: n
    type(mgis_status) :: r
    r = msm_get_uniform_material_property(&
         v, s%ptr, convert_fortran_string(n))
    if( r % exit_status .eq. MGIS_SUCCESS) then
       v = ieee_value(v, ieee_quiet_nan)
    end if
  end function material_state_manager_get_uniform_material_property
  !
  function material_state_manager_get_non_uniform_material_property(v, s, n)&
       result(r)
    use, intrinsic :: iso_c_binding, only: c_f_pointer
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function msm_get_non_uniform_material_property(v, s, n) &
            bind(c,name = 'mgis_bv_material_state_manager_get_non_uniform_material_property') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_char
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out)      :: v
         type(c_ptr), intent(in),value :: s
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: r
       end function msm_get_non_uniform_material_property
    end interface
    real(kind=8), dimension(:), pointer, intent(out) :: v
    type(MaterialStateManager), intent(in) :: s
    character(len=*), intent(in) :: n
    type(mgis_status) :: r
    type(c_ptr) values
    integer nig
    nullify(v)
    r = material_state_manager_get_number_of_integration_points(nig, s)
    if( r%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    r = msm_get_non_uniform_material_property(&
         values, s%ptr, convert_fortran_string(n))
    if( r % exit_status .eq. MGIS_SUCCESS) then
       call c_f_pointer(values, v, [nig])
    end if
  end function material_state_manager_get_non_uniform_material_property
  !
  function material_state_manager_get_internal_state_variables_stride(n_isvs, s) &
       result(r)
    use, intrinsic :: iso_c_binding, only: c_size_t
    use mgis, only: mgis_status, MGIS_SUCCESS
    use mgis_fortran_utilities
    implicit none
    interface
       function msm_get_internal_state_variables_stride_wrapper(nc_isvs, s) &
            bind(c,name = 'mgis_bv_material_state_manager_get_internal_state_variables_stride') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         integer(kind=c_size_t), intent(out) :: nc_isvs
         type(c_ptr), intent(in),value :: s
         type(mgis_status) :: r
       end function msm_get_internal_state_variables_stride_wrapper
    end interface
    integer, intent(out) :: n_isvs
    type(MaterialStateManager), intent(in) :: s
    type(mgis_status) :: r
    integer(kind=c_size_t) nc_isvs
    r = msm_get_internal_state_variables_stride_wrapper(nc_isvs, s %ptr)
    if( r%exit_status .eq. MGIS_SUCCESS) then
       n_isvs = nc_isvs
    else
       n_isvs = -1
    end if
  end function material_state_manager_get_internal_state_variables_stride
  !
  function material_state_manager_get_internal_state_variables(isvs, s) &
       result(r)
    use, intrinsic :: iso_c_binding, only: c_size_t, c_f_pointer
    use mgis, only: mgis_status, MGIS_SUCCESS
    use mgis_fortran_utilities
    implicit none
    interface
       function msm_get_internal_state_variables_wrapper(isvs_ptr, s) &
            bind(c,name = 'mgis_bv_material_state_manager_get_internal_state_variables') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         use mgis, only: mgis_status
         type(c_ptr), intent(out) :: isvs_ptr
         type(c_ptr), intent(in),value :: s
         type(mgis_status) :: r
       end function msm_get_internal_state_variables_wrapper
    end interface
    real(kind=8), dimension(:,:), pointer, intent(out) :: isvs
    type(MaterialStateManager), intent(in) :: s
    type(mgis_status) :: r
    type(c_ptr) isvs_ptr
    integer n_isvs
    integer n
    nullify(isvs)
    r = msm_get_internal_state_variables_wrapper(isvs_ptr, s %ptr)
    if( r%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    r = material_state_manager_get_internal_state_variables_stride(n_isvs, s)
    if( r%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    r = material_state_manager_get_number_of_integration_points(n, s)
    if( r%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    call c_f_pointer(isvs_ptr, isvs, [n_isvs, n])
  end function material_state_manager_get_internal_state_variables
  !
  function material_state_manager_set_uniform_external_state_variable(s, n, v) &
       result(r)
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    implicit none
    interface
       function msm_set_uniform_external_state_variable_wrapper(s, n, v) &
            bind(c,name = 'mgis_bv_material_state_manager_set_uniform_external_state_variable') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_char
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(in),value :: s
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         real(kind=c_double), intent(in), value :: v
         type(mgis_status) :: r
       end function msm_set_uniform_external_state_variable_wrapper
    end interface
    type(MaterialStateManager), intent(in) :: s
    character(len=*), intent(in) :: n
    real(kind=8),     intent(in) :: v
    type(mgis_status) :: r
    r = msm_set_uniform_external_state_variable_wrapper(s %ptr, &
         convert_fortran_string(n), v)
  end function material_state_manager_set_uniform_external_state_variable
  !
  function material_state_manager_is_external_state_variable_defined(b, s, n)&
       result(r)
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function msm_is_external_state_variable_defined(b, s, n) &
            bind(c,name = 'mgis_bv_material_state_manager_is_external_state_variable_defined') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_char
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_int), intent(out) :: b
         type(c_ptr), intent(in),value :: s
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: r
       end function msm_is_external_state_variable_defined
    end interface
    logical, intent(out) :: b
    type(MaterialStateManager), intent(in) :: s
    character(len=*), intent(in) :: n
    type(mgis_status) :: r
    integer bc
    r = msm_is_external_state_variable_defined(&
         bc, s%ptr, convert_fortran_string(n))
    if( r % exit_status .eq. MGIS_SUCCESS) then
       b = bc .eq. 1
    else
       b = .false.
    end if
  end function material_state_manager_is_external_state_variable_defined
  !
  function material_state_manager_is_external_state_variable_uniform(b, s, n)&
       result(r)
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function msm_is_external_state_variable_uniform(b, s, n) &
            bind(c,name = 'mgis_bv_material_state_manager_is_external_state_variable_uniform') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_char
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_int), intent(out) :: b
         type(c_ptr), intent(in),value :: s
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: r
       end function msm_is_external_state_variable_uniform
    end interface
    logical, intent(out) :: b
    type(MaterialStateManager), intent(in) :: s
    character(len=*), intent(in) :: n
    type(mgis_status) :: r
    integer bc
    r = msm_is_external_state_variable_uniform(&
         bc, s%ptr, convert_fortran_string(n))
    if( r % exit_status .eq. MGIS_SUCCESS) then
       b = bc .eq. 1
    else
       b = .false.
    end if
  end function material_state_manager_is_external_state_variable_uniform
  !
  function material_state_manager_get_uniform_external_state_variable(v, s, n)&
       result(r)
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function msm_get_uniform_external_state_variable(v, s, n) &
            bind(c,name = 'mgis_bv_material_state_manager_get_uniform_external_state_variable') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_char
         use mgis, only: mgis_status
         implicit none
         real(kind=c_double), intent(out) :: v
         type(c_ptr), intent(in),value :: s
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: r
       end function msm_get_uniform_external_state_variable
    end interface
    real(kind=8), intent(out) :: v
    type(MaterialStateManager), intent(in) :: s
    character(len=*), intent(in) :: n
    type(mgis_status) :: r
    r = msm_get_uniform_external_state_variable(&
         v, s%ptr, convert_fortran_string(n))
    if( r % exit_status .eq. MGIS_SUCCESS) then
       v = ieee_value(v, ieee_quiet_nan)
    end if
  end function material_state_manager_get_uniform_external_state_variable
  !
  function material_state_manager_get_non_uniform_external_state_variable(v, s, n)&
       result(r)
    use, intrinsic :: iso_c_binding, only: c_f_pointer
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use mgis_fortran_utilities
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function msm_get_non_uniform_external_state_variable(v, s, n) &
            bind(c,name = 'mgis_bv_material_state_manager_get_non_uniform_external_state_variable') &
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_char
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out)      :: v
         type(c_ptr), intent(in),value :: s
         character(len=1,kind=c_char), dimension(*), intent(in) :: n
         type(mgis_status) :: r
       end function msm_get_non_uniform_external_state_variable
    end interface
    real(kind=8), dimension(:), pointer, intent(out) :: v
    type(MaterialStateManager), intent(in) :: s
    character(len=*), intent(in) :: n
    type(mgis_status) :: r
    type(c_ptr) values
    integer nig
    nullify(v)
    r = material_state_manager_get_number_of_integration_points(nig, s)
    if( r%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    r = msm_get_non_uniform_external_state_variable(&
         values, s%ptr, convert_fortran_string(n))
    if( r % exit_status .eq. MGIS_SUCCESS) then
       call c_f_pointer(values, v, [nig])
    end if
  end function material_state_manager_get_non_uniform_external_state_variable
  ! \brief return the transpose of the tangent operator
  function material_data_manager_get_tangent_operator(K, m) result(r)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
    use mgis, only: mgis_status, MGIS_SUCCESS
    implicit none
    interface
       function material_data_manager_get_tangent_operator_wrapper(K, m) &
            bind(c,name = 'mgis_bv_material_data_manager_get_tangent_operator')&
            result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: K
         type(c_ptr), intent(in), value :: m
         type(mgis_status) :: r
       end function material_data_manager_get_tangent_operator_wrapper
    end interface
    real(kind=8), dimension(:,:,:), pointer, intent(out) :: K
    type(MaterialDataManager), intent(in) :: m
    type(mgis_status) :: r
    type(MaterialStateManager) :: s0
    type(c_ptr) :: p
    integer n
    integer gs
    integer ths
    nullify(K)
    r = material_data_manager_get_state_0(s0, m)
    if( r%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    r = material_data_manager_get_tangent_operator_wrapper(p, m %ptr)
    if( r%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    r = material_state_manager_get_gradients_stride(gs, s0)
    if( r%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    r = material_state_manager_get_thermodynamic_forces_stride(gs, s0)
    if( r%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    r = material_state_manager_get_number_of_integration_points(n, s0)
    if( r%exit_status .ne. MGIS_SUCCESS) then
       return
    end if
    if( r%exit_status .eq. MGIS_SUCCESS) then
       call c_f_pointer(p,K,[ths,gs, n])
    end if
  end function material_data_manager_get_tangent_operator
  !
  function free_material_data_manager(ptr) result(r)
    use, intrinsic :: iso_c_binding, only: c_associated
    use mgis
    implicit none
    interface
       function free_material_data_manager_wrapper(ptr) &
            bind(c, name='mgis_bv_free_material_data_manager') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis
         implicit none
         type(c_ptr), intent(inout) :: ptr
         type(mgis_status) :: r
       end function free_material_data_manager_wrapper
    end interface
    type(MaterialDataManager), intent(inout) :: ptr
    type(mgis_status) :: r
    if (c_associated(ptr%ptr)) then
       r = free_material_data_manager_wrapper(ptr%ptr)
    end if
  end function free_material_data_manager
  !
  function update_material_data_manager(d) result(r)
    use mgis
    implicit none
    interface
       function update_material_data_manager_wrapper(ptr) &
            bind(c, name='mgis_bv_update_material_data_manager') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis
         implicit none
         type(c_ptr), intent(in), value :: ptr
         type(mgis_status) :: r
       end function update_material_data_manager_wrapper
    end interface
    type(MaterialDataManager), intent(in) :: d
    type(mgis_status) :: r
    r = update_material_data_manager_wrapper(d%ptr)
  end function update_material_data_manager
  !
  function revert_material_data_manager(d) result(r)
    use mgis
    implicit none
    interface
       function revert_material_data_manager_wrapper(ptr) &
            bind(c, name='mgis_bv_revert_material_data_manager') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis
         implicit none
         type(c_ptr), intent(in), value :: ptr
         type(mgis_status) :: r
       end function revert_material_data_manager_wrapper
    end interface
    type(MaterialDataManager), intent(in) :: d
    type(mgis_status) :: r
    r = revert_material_data_manager_wrapper(d%ptr)
  end function revert_material_data_manager
  ! integrate
  function integrate(r, d, b) result(s)
    use mgis, only: mgis_status
    implicit none
    interface
       function integrate_wrapper(r, d, b) &
            bind(c,name = 'mgis_bv_integrate_2') result(s)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_int
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_int), intent(out) :: r
         type(c_ptr), intent(in),value :: d
         type(c_ptr), intent(in),value :: b
         type(mgis_status) :: s
       end function integrate_wrapper
    end interface
    integer, intent(out) :: r
    type(BehaviourData), intent(in) :: d
    type(Behaviour),     intent(in) :: b
    type(mgis_status) :: s
    s = integrate_wrapper(r, d%ptr, b%ptr)
  end function integrate
  ! 
  function integrate_material_data_manager(r, p, m, i, dt) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t
    use mgis_fortran_utilities, only: convert_to_c_index
    use mgis, only: ThreadPool, mgis_status, report_failure
    implicit none
    interface
       function integrate_material_data_manager_wrapper(r, p, m, i, dt) &
            bind(c,name = 'mgis_bv_integrate_material_data_manager') &
            result(s)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_int), intent(out) :: r
         type(c_ptr), intent(in),value :: p
         type(c_ptr), intent(in),value :: m
         integer,     intent(in),value :: i
         real(kind = c_double), intent(in),value :: dt
         type(mgis_status) :: s
       end function integrate_material_data_manager_wrapper
    end interface
    integer, intent(out) :: r
    type(ThreadPool),          intent(in) :: p
    type(MaterialDataManager), intent(in) :: m
    integer,                   intent(in) :: i
    real(kind = 8),            intent(in) :: dt
    type(mgis_status) :: s
    s = integrate_material_data_manager_wrapper(r, p%ptr, m%ptr, i, dt)
  end function integrate_material_data_manager
  !
  function integrate_material_data_manager_part(r, m, i, dt, ni, ne) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t
    use mgis_fortran_utilities, only: convert_to_c_index
    use mgis, only: mgis_status, report_failure
    implicit none
    interface
       function integrate_material_data_manager_part_wrapper(r, m, i, dt, ni, ne) &
            bind(c,name = 'mgis_bv_integrate_material_data_manager_part') &
            result(s)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t, c_int, c_double
         use mgis, only: mgis_status
         implicit none
         integer(kind=c_int), intent(out) :: r
         type(c_ptr), intent(in),value :: m
         integer,     intent(in),value :: i
         real(kind = c_double), intent(in),value :: dt
         integer(kind = c_size_t), intent(in),value :: ni
         integer(kind = c_size_t), intent(in),value :: ne
         type(mgis_status) :: s
       end function integrate_material_data_manager_part_wrapper
    end interface
    integer, intent(out) :: r
    type(MaterialDataManager), intent(in) :: m
    integer,     intent(in),value :: i
    real(kind = 8),            intent(in) :: dt
    integer :: ni
    integer :: ne
    type(mgis_status) :: s
    integer(kind=c_size_t) :: nic
    integer(kind=c_size_t) :: nec
    if(.not. convert_to_c_index(nic, ni)) then
       s = report_failure("invalid index")
       return
    end if
    if(.not. convert_to_c_index(nec, ne)) then
       s = report_failure("invalid index")
       return
    end if
    nec = nec + 1
    s = integrate_material_data_manager_part_wrapper(r, m%ptr, i, dt, nic, nec)
  end function integrate_material_data_manager_part
end module  mgis_behaviour
