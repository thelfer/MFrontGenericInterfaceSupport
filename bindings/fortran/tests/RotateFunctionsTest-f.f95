subroutine test()
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  type(Behaviour) b
  real(kind=8), dimension(6) :: g ! gradient
  real(kind=8), dimension(9) :: r ! rotation matrix
  logical :: res
  ! reference values of the equivalent plastic strain
  g = (/ 1, 0, 0, 0, 0, 0/)
  r = (/ 0, 1, 0, &
         1, 0, 0, &
         0, 0, 1/)
  call check_status(load_behaviour(b, &
       get_mfront_behaviour_test_library_path(), &
       'OrthotropicElasticity', 'Tridimensional'))
  call check_status(rotate_gradients_in_place(g, b, r))
  res = check(abs(g(1)-0) < 1.d-12, &
       'invalid value of the gradient')
  res = check(abs(g(2)-1) < 1.d-12, &
       'invalid value of the gradient')
  res = check(abs(g(3)-0) < 1.d-12, &
       'invalid value of the gradient')
  res = check(abs(g(4)-0) < 1.d-12, &
       'invalid value of the gradient')
  res = check(abs(g(5)-0) < 1.d-12, &
       'invalid value of the gradient')
  res = check(abs(g(6)-0) < 1.d-12, &
       'invalid value of the gradient')
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
