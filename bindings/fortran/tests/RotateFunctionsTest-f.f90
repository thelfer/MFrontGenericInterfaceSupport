subroutine test()
  use, intrinsic :: iso_c_binding, only: c_size_t
  use mgis
  use mgis_behaviour
  use mgis_testing_utilities
  implicit none
  type(Behaviour) b
  real(kind=8), dimension(9) :: r ! rotation matrix
  real(kind=8), dimension(4) :: ge, me ! gradients in the global and material frames
  real(kind=8), dimension(4) :: v ! rotated values
  real(kind=8), dimension(8) :: v2 ! values at two integration points
  real(kind=8), dimension(4) :: sm, sg ! thermodynamic forces
  real(kind=8), dimension(16) :: Km ! tangent operator in the material frame
  real(kind=8), dimension(32) :: K2 ! tangent operators at two integration points
  integer(kind=c_size_t), parameter :: n = 2 ! number of integration points
  real(kind=8) :: c, s, cste, exy
  integer :: i, j, p
  logical :: res
  ! rotation of 30 degrees around z: r stores row by row the matrix
  ! R = {{c, -s, 0}, {s, c, 0}, {0, 0, 1}}, which is not symmetric (see
  ! RotateFunctionsTest.cxx)
  c = cos(acos(-1.d0) / 6)
  s = sin(acos(-1.d0) / 6)
  cste = sqrt(2.d0)
  r = (/ c, -s, 0.d0, &
         s, c, 0.d0, &
         0.d0, 0.d0, 1.d0/)
  call check_status(load_behaviour(b, &
       get_mfront_behaviour_test_library_path(), &
       'OrthotropicElasticity', 'GeneralisedPlaneStrain'))
  ! gradient with a shear component (stored as sqrt(2) times exy)
  ge = (/ 1.d0, 2.d0, 3.d0, 4.d0 /)
  me = ge
  call check_status(rotate_gradients_in_place(me, b, r))
  ! reference: R^T E R, the gradient in the material frame
  exy = ge(4) / cste
  res = check(abs(me(1) - (c * c * ge(1) + 2 * c * s * exy + s * s * ge(2))) < 1.d-12, &
       'invalid value of the gradient')
  res = check(abs(me(2) - (s * s * ge(1) - 2 * c * s * exy + c * c * ge(2))) < 1.d-12, &
       'invalid value of the gradient')
  res = check(abs(me(3) - ge(3)) < 1.d-12, &
       'invalid value of the gradient')
  res = check(abs(me(4) - cste * (c * s * (ge(2) - ge(1)) + (c * c - s * s) * exy)) < 1.d-12, &
       'invalid value of the gradient')
  ! the other rotation functions of the gradients must give the same results
  call check_status(rotate_gradients_out_of_place(v, b, ge, r))
  res = check(all(abs(v - me) < 1.d-12), 'invalid value of the gradient')
  ! the rotation of the thermodynamic forces is the inverse rotation
  call check_status(rotate_array_of_thermodynamic_forces_out_of_place(v2, b, (/ me, me /), r, n))
  res = check(all(abs(v2 - (/ ge, ge /)) < 1.d-12), 'invalid value of the thermodynamic forces')
  ! tangent operator, stored row by row, with no symmetry: the rotated
  ! operator Kg must satisfy Kg ge = rotate_thermodynamic_forces(Km me)
  do i = 1, 16
     Km(i) = real(i, 8) + 0.1d0 * real(mod(7 * i, 11), 8)
  end do
  do i = 1, 4
     sm(i) = 0
     do j = 1, 4
        sm(i) = sm(i) + Km((i - 1) * 4 + j) * me(j)
     end do
  end do
  call check_status(rotate_thermodynamic_forces_out_of_place(sg, b, sm, r))
  K2 = (/ Km, Km /)
  call check_status(rotate_array_of_tangent_operator_blocks_in_place(K2, b, r, n))
  do p = 0, 1
     do i = 1, 4
        v(i) = 0
        do j = 1, 4
           v(i) = v(i) + K2(p * 16 + (i - 1) * 4 + j) * ge(j)
        end do
     end do
     res = check(all(abs(v - sg) < 1.d-12 * max(1.d0, maxval(abs(sg)))), &
          'invalid value of the tangent operators')
  end do
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
