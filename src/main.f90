program second_order
  use settings
  use solver_utils
  use dim_flexible_multi_body_dyn
  implicit none
  
  integer, parameter        :: t0=0, t_final=1
  integer, parameter        :: nSteps=10 ! number of time steps
  integer, parameter        :: nVars=1
  integer                   :: i,j,k
  real(dp), parameter       :: h = dble(t_final-t0)/dble(nSteps) ! time step size
  real(dp)                  :: q(0:nSteps,nVars)
  real(dp)                  :: q_dot(0:nSteps,1:nVars)
  real(dp)                  :: q_double_dot(0:nSteps, 1:nVars)  ! state variable, its first and second derivatives

  write(*,*) "-----------------------"
  write(*,*) "Second order ODE solver"
  write(*,*) "-----------------------"

  !-------------------
  ! Initial conditions
  !-------------------
  q(0,:) = (/ 1.0_dp /)           ! set initial condition for q(x=1.0) 
  q_dot(0,:) =  (/ 1.0_dp /)      ! set initial condition for q_dot

  !reshape((/ (i, i = 1,nvars**nvars) /),shape(q_dot(0,:,:)))    ! set initial condition for q_dot

  print *,h,  q(0,:), q_dot(0,:)
  print *,  get_extrapolated_q(h,  q(0,:), q_dot(0,:))
stop
  !  print*,q_dot(0,:,:)
  !  write(*,*) "The initial conditions are:"
  !  write(*,*) "q     =", q(0,:)
  !  write(*,*) "q_dot =", (  q_dot(0,:,i), i= 1,nVars )
  call test_get_extrapolated_q
!  call test04
  print*,  get_bdf_coeffs(3,2)

end program second_order

! take a simple multi body equation and use it for implementation
! fetch routine to get coeffs for the approximation of first and second derivatives
! swap out solvers if needed
! treat q as scalar if needed

!----------------------------------------------------------------------
! routine to test the extrapolation of the first and second derivatives
!----------------------------------------------------------------------
subroutine test_get_extrapolated_q
use settings
use solver_utils
implicit none

real(dp):: q(2)
real(dp) ::q0(2),q0_dot(2), q0_double_dot(2)
real(dp):: x ,delta_t

! set constants
delta_t = 0.01_dp
x=1._dp

! give starting values for extrapolation
q0=(/ x**3, sin(x) /)
q0_dot = (/ 3._dp*x**2, cos(x) /)
q0_double_dot = (/ 6._dp*x, -sin(x) /)

! call the extrapolation routine
print*, "old q=", q0
q(:) = get_extrapolated_q(delta_t,  q0, q0_dot, q0_double_dot)
print*, "new q=", q

! call the extrapolation routine
print*, "old q=", q0
q(:) = get_extrapolated_q(delta_t, q0, q0_dot)
print*, "new q=", q

return
end subroutine test_get_extrapolated_q
