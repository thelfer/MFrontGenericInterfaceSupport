module mgis
  use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_char, c_null_ptr
  use, intrinsic :: iso_fortran_env, only: int64
  use mgis_fortran_utilities
  ! enumeration
  enum , bind(c) 
     enumerator :: MGIS_SUCCESS = 0
     enumerator :: MGIS_FAILURE = 1
  end enum
  type, bind(c) :: mgis_status
     integer :: exit_status
     type(c_ptr) :: msg = c_null_ptr
  end type mgis_status
  type :: ThreadPool
     type(c_ptr) :: ptr
  end type ThreadPool
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
  !
  function create_thread_pool(p, n) result(s)
    use, intrinsic :: iso_c_binding, only: c_size_t, c_associated
    use mgis_fortran_utilities
    implicit none
    interface
       function create_thread_pool_wrapper(p, n) &
            bind(c,name = 'mgis_create_thread_pool') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
         import mgis_status
         implicit none
         type(c_ptr), intent(out) :: p
         integer(kind=c_size_t), intent(in), value :: n
         type(mgis_status) :: r
       end function create_thread_pool_wrapper
    end interface
    type(ThreadPool), intent(out) :: p
    integer, intent(in) :: n
    type(mgis_status) :: s
    integer(kind=c_size_t) nc
    if (n.lt.1) then
       s = report_failure('invalid number of threads')
       return
    end if
    nc = n
    s = create_thread_pool_wrapper(p%ptr, nc)
  end function create_thread_pool
  !
  function free_thread_pool(p) result(r)
    use, intrinsic :: iso_c_binding, only: c_associated
    implicit none
    interface
       function free_thread_pool_wrapper(p) &
            bind(c, name='mgis_free_thread_pool') result(r)
         use, intrinsic :: iso_c_binding, only: c_ptr
         import mgis_status
         implicit none
         type(c_ptr), intent(inout) :: p
         type(mgis_status) :: r
       end function free_thread_pool_wrapper
    end interface
    type(ThreadPool), intent(inout) :: p
    type(mgis_status) :: r
    if (c_associated(p%ptr)) then
       r = free_thread_pool_wrapper(p%ptr)
    end if
  end function free_thread_pool
end module mgis
