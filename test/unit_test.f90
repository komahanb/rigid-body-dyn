

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
