program second_order
  use settings
  use dim_flexible_multi_body_dyn
  implicit none

  integer, parameter        :: t0=0, t_final=1
  integer, parameter        :: nSteps=10 ! number of time steps
  integer, parameter        :: nVars=3
  integer                   :: i,j,k
  real(dp), parameter       :: h = dble(t_final-t0)/dble(nSteps) ! time step size
  real(dp)                  :: q(0:nSteps,nVars)
  real(dp)                  :: q_dot(0:nSteps,1:nVars,1:nVars)
  real(dp)                  :: q_double_dot(0:nSteps, 1:nVars,1:nVars,1:nVars)  ! state variable, its first and second derivatives

  write(*,*) "-----------------------"
  write(*,*) "Second order ODE solver"
  write(*,*) "-----------------------"

  !-------------------
  ! Initial conditions
  !-------------------
  q(0,:) = (/ (i, i = 1, nVars) /)                        ! set initial condition for q
  q_dot(0,:,:) = reshape((/ (i, i = 1,nvars**nvars) /),shape(q_dot(0,:,:)))    ! set initial condition for q_dot
  !  print*,q_dot(0,:,:)
  !  write(*,*) "The initial conditions are:"
  !  write(*,*) "q     =", q(0,:)
  !  write(*,*) "q_dot =", (  q_dot(0,:,i), i= 1,nVars )
  call test_get_extrapolated_q

end program second_order

! take a simple multi body equation and use it for implementation
! fetch routine to get coeffs for the approximation of first and second derivatives
! swap out solvers if needed
! treat q as scalar if needed

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
q(:) = get_extrapolated_q(dT=delta_t, old_q = q0, old_q_dot=q0_dot, old_q_double_dot=q0_double_dot)
print*, "new q=", q

return
end subroutine test_get_extrapolated_q
