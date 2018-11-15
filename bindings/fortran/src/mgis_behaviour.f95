module mgis_behaviour
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
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
  type :: behaviour
    private
    type(c_ptr) :: ptr = c_null_ptr
  end type behaviour
contains
  function load_behaviour(b,l,bn,h) result(s)
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    implicit none
    interface
       function load_behaviour_wrapper(ptr,l,b,h) bind(c,name = 'mgis_bv_load_behaviour') result(r)
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
  ! behaviour_get_library
  function behaviour_get_library(l,b) result(s)
    use mgis_fortran_utilities
    use mgis, only: mgis_status
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
    l = convert_c_string(lp)
  end function behaviour_get_library
  ! behaviour_get_behaviour_name
  function behaviour_get_behaviour_name(l,b) result(s)
    use mgis_fortran_utilities
    use mgis, only: mgis_status
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
    l = convert_c_string(lp)
  end function behaviour_get_behaviour_name
  ! behaviour_get_source
  function behaviour_get_source(l,b) result(s)
    use mgis_fortran_utilities
    use mgis, only: mgis_status
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
    l = convert_c_string(lp)
  end function behaviour_get_source
  ! behaviour_get_hypothesis
  function behaviour_get_hypothesis(l,b) result(s)
    use mgis_fortran_utilities
    use mgis, only: mgis_status
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
    l = convert_c_string(lp)
  end function behaviour_get_hypothesis
  ! behaviour_get_tfel_version
  function behaviour_get_tfel_version(l,b) result(s)
    use mgis_fortran_utilities
    use mgis, only: mgis_status
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
    l = convert_c_string(lp)
  end function behaviour_get_tfel_version
  ! free behaviour
  ! Use a finalizer to do automatic cleanup off the C structures.
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
    type(behaviour), intent(inout) :: ptr
    type(mgis_status) :: r
    if (c_associated(ptr%ptr)) then
       r = free_behaviour_wrapper(ptr%ptr)
    end if
  end function free_behaviour
end module  mgis_behaviour
