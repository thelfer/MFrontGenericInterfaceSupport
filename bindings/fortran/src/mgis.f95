module mgis
  use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_char
  use, intrinsic :: iso_fortran_env, only: int64
  use mgis_fortran_utilities
  ! enumeration
  enum , bind(c) 
     enumerator :: MGIS_SUCCESS = 0
     enumerator :: MGIS_FAILURE = 1
  end enum
  type, bind(c) :: mgis_status
     integer :: exit_status
     type(c_ptr) :: msg
  end type mgis_status
  interface
     function report_success() bind(c,name = 'mgis_report_success') result(r)
       import mgis_status
       implicit none
       type(mgis_status) :: r
     end function report_success
  end interface
contains
  function report_failure(msg) result(r)
    use, intrinsic :: iso_c_binding, only: c_ptr, c_char,c_loc
    implicit none
    interface
       function report_failure_wrapper(msg) bind(c,name = 'mgis_report_failure') result(r)
         import mgis_status, c_ptr, c_char
         implicit none
         character(len=1,kind=c_char), dimension(*), intent(in) :: msg
         type(mgis_status) :: r
       end function report_failure_wrapper
    end interface
    character(len=*), intent(in) :: msg
    type(mgis_status) :: r
    character(len=1,kind=c_char) :: s(len(msg)+1)
    s = convert_fortran_string(msg)
    r = report_failure_wrapper(s)
  end function report_failure
  function get_error_message(s) result(msg)
    implicit none
    type(mgis_status) :: s
    character(len=:), allocatable :: msg
    msg = convert_c_string(s%msg)
  end function get_error_message
end module mgis
