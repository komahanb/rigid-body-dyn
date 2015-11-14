! program to test scalar_class
program test_scalar

  ! program references
  use scalar_class, only : scalar
  use global_constants, only : dp

  ! program settings
  implicit none

  ! program variables
  class(scalar), pointer :: z
  type(scalar)  :: z2
  real(dp)      :: x, y

  z2 = scalar(x=1.5_dp)
  call z2%print()

  call z2%set_real(1.0_dp)
  call z2%print()

  call z2%set_cplx(1.0_dp)
  call z2%print()

  print*,  z2%get_real()
  print*,  z2%get_cplx()

end program test_scalar

! dummy routines that need impl
subroutine residual
end subroutine residual

subroutine residual2
end subroutine residual2
