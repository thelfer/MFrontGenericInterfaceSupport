module mgis_testing_utilities
  logical :: status = .true.
  integer :: number_of_tests = 0
  integer :: number_of_successes = 0
  integer :: number_of_failures = 0
contains
  !
  function get_mfront_behaviour_test_library_path() result(l)
    use, intrinsic :: ISO_C_BINDING, only: c_ptr
    use mgis_fortran_utilities
    implicit none
    interface
       function mgis_get_mfront_behaviour_test_library_path_wrapper() &
            bind(c,name = 'mgis_get_mfront_behaviour_test_library_path') result(l)
         import c_ptr
         implicit none
         type(c_ptr) :: l
       end function mgis_get_mfront_behaviour_test_library_path_wrapper
    end interface
    character(len=:), allocatable :: l
    l = convert_c_string(mgis_get_mfront_behaviour_test_library_path_wrapper())
  end function get_mfront_behaviour_test_library_path
  !
  function get_mfront_model_test_library_path() result(l)
    use, intrinsic :: ISO_C_BINDING, only: c_ptr
    use mgis_fortran_utilities
    implicit none
    interface
       function mgis_get_mfront_model_test_library_path_wrapper() &
            bind(c,name = 'mgis_get_mfront_model_test_library_path') result(l)
         import c_ptr
         implicit none
         type(c_ptr) :: l
       end function mgis_get_mfront_model_test_library_path_wrapper
    end interface
    character(len=:), allocatable :: l
    l = convert_c_string(mgis_get_mfront_model_test_library_path_wrapper())
  end function get_mfront_model_test_library_path
  !
  function get_tfel_version() result(l)
    use, intrinsic :: ISO_C_BINDING, only: c_ptr
    use mgis_fortran_utilities
    implicit none
    interface
       function mgis_get_tfel_version_wrapper() &
            bind(c,name = 'mgis_get_tfel_version') result(l)
         import c_ptr
         implicit none
         type(c_ptr) :: l
       end function mgis_get_tfel_version_wrapper
    end interface
    character(len=:), allocatable :: l
    l = convert_c_string(mgis_get_tfel_version_wrapper())
  end function get_tfel_version
  ! 
  function check (c, msg)
    implicit none
    logical, intent(in) :: c
    character(len=*), intent(in) :: msg
    logical :: check
    if(.not. c) then
       write(*,*) msg
       status = .false.
       number_of_failures = number_of_failures + 1
    else
       number_of_successes = number_of_successes + 1
    end if
    number_of_tests = number_of_tests + 1
    check = c
  end function check
  !
  function check_string (s1, s2, msg) result(c)
    implicit none
    character(len=*), intent(in) :: s1
    character(len=*), intent(in) :: s2
    character(len=*), intent(in) :: msg
    logical :: c
    c = check(s1 .eq. s2, msg)
  end function check_string
  !
  subroutine check_status(s)
    use mgis, only: mgis_status, MGIS_SUCCESS, get_error_message
    implicit none
    type(mgis_status), intent(in) :: s
    if ( s%exit_status .ne. MGIS_SUCCESS) then
       write(*,*) "invalid function call: ", get_error_message(s)
       stop -1
    end if
  end subroutine check_status
  !
  subroutine tests_summary()
    implicit none
    write(*,*) '---------------------------------------------------'
    if(status) then
       write (*,*) 'SUCCESS'
    else
       write (*,*) 'FAILURE'
    endif
    write (*,*) "Total number of tests:",  number_of_tests
    if(.not. status) then
       write (*,*) "Number of succeeded tests:",  number_of_successes, &
            '(',(number_of_failures*100.)/number_of_tests,'%)'
       write (*,*) "Number of failed tests:",  number_of_failures, &
            '(',(number_of_failures*100.)/number_of_tests,'%)'
    end if
    write(*,*) '---------------------------------------------------'
  end subroutine tests_summary
end module mgis_testing_utilities
