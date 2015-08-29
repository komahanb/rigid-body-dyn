program second_order
  use constants
  use solver_utils
  use dim_flexible_multi_body_dyn
  implicit none

  integer, parameter        :: t0=0, t_final=1
  integer, parameter        :: nSteps=10 ! number of time steps
  integer                   :: i,j,k
  real(dp)                  :: q(0:nSteps,dim_q)
  real(dp)                  :: q_dot(0:nSteps,dim_q)
  real(dp)                  :: q_double_dot(0:nSteps, dim_q)  ! state variable, its first and second derivatives

  !-------------------
  ! Initial conditions
  !-------------------

  q(0,:)     = (/ 1.0_dp, 1.0_dp /)           ! set initial condition for q(x=1.0) 
  q_dot(0,:) =  (/ 1.0_dp, 1.0_dp /)      ! set initial condition for q_dot

  call test_get_extrapolated_q


end program second_order

!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine test_
  use constants
  use solver_utils
  implicit none

end subroutine test_

!----------------------------------------------------------------------

!----------------------------------------------------------------------
subroutine test
  use constants
  use solver_utils
  implicit none

end subroutine test

!----------------------------------------------------------------------
! routine to test the extrapolation of the first and second derivatives
!----------------------------------------------------------------------
subroutine test_get_extrapolated_q
  use constants
  use solver_utils
  implicit none

  real(dp):: q(2)
  real(dp) ::q0(2),q0_dot(2), q0_double_dot(2)
  real(dp):: x

  print*, "------------------------------------------------------"
  print*, "---------test_get_extrapolated_q----------------------"
  print*, "------------------------------------------------------"

  ! set constants
  x=1._dp

  ! give starting values for extrapolation
  q0=(/ x**3, sin(x) /)
  q0_dot = (/ 3._dp*x**2, cos(x) /)
  q0_double_dot = (/ 6._dp*x, -sin(x) /)

  ! call the extrapolation routine
  print*, "old q=", q0
  q(:) = get_extrapolated_q( q0, q0_dot, q0_double_dot)
  print*, "new q=", q

  ! call the extrapolation routine
  print*, "old q=", q0
  q(:) = get_extrapolated_q( q0, q0_dot)
  print*, "new q=", q


  x=1._dp+del_t
  q0=(/ x**3, sin(x) /)

  print*, "ext q=", q0
  print*, "--------------complete--------------------------------"
  print*, ""

  return
end subroutine test_get_extrapolated_q
