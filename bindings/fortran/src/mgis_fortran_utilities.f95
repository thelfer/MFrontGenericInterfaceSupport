module mgis_fortran_utilities
  interface
     ! exposes strlen to fortran
     function strlen_wrapper(s) bind(c,name='strlen') result(r)
       use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
       implicit none
       type(c_ptr),intent(in),value :: s
       integer(kind=c_size_t) :: r
     end function strlen_wrapper
  end interface
contains
  ! convert C string to fortran string
  function convert_c_string (cs) result (s)
    use, intrinsic :: iso_c_binding, only: c_ptr,c_char,c_null_char,&
         c_size_t,c_f_pointer,c_associated
    implicit none
    type(c_ptr), intent(in) :: cs
    character(len=:), allocatable :: s
    character(kind=c_char), dimension(:), pointer :: fptr => null()
    integer(kind=c_size_t) :: i,length
    if(.not. c_associated(cs)) then
       allocate(character(len=0)::s)
    else       
       length = strlen_wrapper(cs)
       allocate(character(len=length)::s)
       call c_f_pointer(cs,fptr,[length])
       do i=1,length
          s(i:i)=fptr(i)
       end do
    end if
  end function convert_c_string
  ! convert fortran string to "c" string
  pure function convert_fortran_string (s) result (cs)
    use, intrinsic :: iso_c_binding, only: c_char, c_null_char
    implicit none
    character(len=*), intent(in) :: s
    character(len=1,kind=c_char) :: cs(len(s)+1)
    integer                      :: n, i
    n = len(s)
    do i = 1, n
       cs(i) = s(i:i)
    end do
    cs(n + 1) = c_null_char
  end function convert_fortran_string
end module mgis_fortran_utilities
