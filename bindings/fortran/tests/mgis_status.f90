subroutine test()
  use mgis
  use mgis_testing_utilities
  implicit none
  type(mgis_status) :: s
  logical :: b
  character(len=30) :: msg = 'failure generated from fortran'
  ! test of successful status
  s = report_success()
  b = check(s%exit_status .eq. MGIS_SUCCESS,'invalid status')
  b = check(get_error_message(s) .eq. '','invalid error message')
  ! test of failed status
  s = report_failure(msg)
  b = check(s%exit_status .eq. MGIS_FAILURE,'invalid status')
  b = check(get_error_message(s) .eq. msg,'invalid error message')
end subroutine test

program main
  use mgis_testing_utilities
  call test()
  call tests_summary()
  if(.not. status) then
     stop -1
  end if
end program main
