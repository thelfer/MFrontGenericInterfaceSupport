module mgis_model
  use mgis_behaviour
  type, extends(Behaviour) :: Model
  end type Model
contains
  !
  function load(m,l,mn,h) result(s)
    use mgis_fortran_utilities
    use mgis, only: mgis_status
    implicit none
    interface
       function load_wrapper(ptr,l,m,h) &
            bind(c,name = 'mgis_model_load') result(r)
         use, intrinsic :: iso_c_binding, only: c_char, c_ptr
         use mgis, only: mgis_status
         implicit none
         type(c_ptr), intent(out) :: ptr
         character(len=1,kind=c_char), dimension(*), intent(in) :: l
         character(len=1,kind=c_char), dimension(*), intent(in) :: m
         character(len=1,kind=c_char), dimension(*), intent(in) :: h
         type(mgis_status) :: r
       end function load_wrapper
    end interface
    type(Model), intent(out) :: m
    character(len=*), intent(in) :: l
    character(len=*), intent(in) :: mn
    character(len=*), intent(in) :: h
    type(mgis_status) :: s
    s = load_wrapper(m%ptr, convert_fortran_string(l), &
         convert_fortran_string(mn), &
         convert_fortran_string(h))
  end function load
  ! free behaviour
  function free_model(ptr) result(r)
    use, intrinsic :: iso_c_binding, only: c_associated
    use mgis
    implicit none
    interface
       function free_model_wrapper(ptr) bind(c, name='mgis_model_free_model') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         use mgis
         implicit none
         type(c_ptr), intent(inout) :: ptr
         type(mgis_status) :: r
       end function free_model_wrapper
    end interface
    type(Model), intent(inout) :: ptr
    type(mgis_status) :: r
    if (c_associated(ptr%ptr)) then
       r = free_model_wrapper(ptr%ptr)
    end if
  end function free_model
end module  mgis_model
