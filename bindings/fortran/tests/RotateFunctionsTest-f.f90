subroutine test()
  use, intrinsic :: iso_c_binding, only: c_size_t
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  type(Behaviour) b
  real(kind=8), dimension(6) :: g ! gradient
  real(kind=8), dimension(9) :: r ! rotation matrix
  real(kind=8), dimension(6) :: g0, v ! unrotated and rotated values
  real(kind=8), dimension(12) :: v2 ! values at two integration points
  real(kind=8), dimension(36) :: K0, Kr ! tangent operators
  real(kind=8), dimension(72) :: K2 ! tangent operators at two integration points
  integer(kind=c_size_t), parameter :: n = 2 ! number of integration points
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
  write(*,*) 'g:', g(1)
  write(*,*) 'g:', g(2)
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
  ! the other rotation functions must give the same results
  g0 = (/ 1, 0, 0, 0, 0, 0/)
  v = g0
  call check_status(rotate_gradients_out_of_place(v, b, g0, r))
  res = check(all(abs(v-g) < 1.d-12), 'invalid value of the gradient')
  v2 = (/ g0, g0 /)
  call check_status(rotate_array_of_thermodynamic_forces_out_of_place(v2, b, (/ g0, g0 /), r, n))
  res = check(all(abs(v2-(/ g, g /)) < 1.d-12), 'invalid value of the thermodynamic forces')
  ! the rotation exchanges the first two axes: K0(1,1) is moved to Kr(2,2)
  K0 = 0
  K0(1) = 1
  Kr = 0
  Kr(8) = 1
  K2 = (/ K0, K0 /)
  call check_status(rotate_array_of_tangent_operator_blocks_in_place(K2, b, r, n))
  res = check(all(abs(K2-(/ Kr, Kr /)) < 1.d-12), 'invalid value of the tangent operators')
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
