!*********************************************************************|
! Example calling program for the finite_diff module
!*********************************************************************|
program test_finite_diff
  use finite_diff,only:finite_difference
  implicit none
  integer, parameter                    :: sp = kind(0.0) 
  integer(sp), parameter                :: dp = kind(0.0d0)
  real(dp)::q(2),jac_approx(2,2), jac_exact(2,2)

  ! get a random input variable
  call random_number(q)
  
  ! find the approximate jacobian
  jac_approx = finite_difference(q,1.0d-6)

  ! find the actual jacobian
  call exact_diff_test1(q,jac_exact)  
  
  ! find the max difference
  write(*,*) "Max Diff : ", maxval(jac_approx - jac_exact)

end program test_finite_diff

!*********************************************************************!
! Exact JAC for TEST 1
!*********************************************************************!
subroutine exact_diff_test1(q, jac)
  implicit none

  integer, parameter                    :: sp = kind(0.0)  
  integer(sp), parameter                :: dp = kind(0.0d0)
  real(dp) :: q(2), jac(2,2)

  jac(1,1) = exp(q(2))*cos(q(1)) 
  jac(1,2) = sin(q(1))*exp(q(2))
  jac(2,1) = 0.0_dp
  jac(2,2) = -sin(q(2))

end subroutine exact_diff_test1

!*********************************************************************|
! Example implementation of the function call to compute FD derivative
! This function must be implemented by the calling program when FD 
! jacobian is wanted. 
!*********************************************************************|
subroutine residual(q, f)

  implicit none

  integer, parameter                    :: sp = kind(0.0)  
  integer(sp), parameter                :: dp = kind(0.0d0)

  real(dp), intent(in), dimension(2)  :: q ! the variables
  real(dp), intent(out), dimension(2) :: f ! the vector valued func

  f = (/ sin(q(1))*exp(q(2)), cos(q(2)) /)

end subroutine residual
