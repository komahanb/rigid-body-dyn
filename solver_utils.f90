module solver_utils
use settings
implicit none

private
public get_extrapolated_q

contains

!--------------------------------------------------------------------------
! returns the extrapolated value of q based on first and second derivatives
!-------------------------------------------------------------------------
function get_extrapolated_q(dT, old_q, old_q_dot, old_q_double_dot) result(new_q)
real(dp), intent(in)           :: dT
real(dp), intent(in)           :: old_q(:), old_q_dot(:)
real(dp), intent(in), optional :: old_q_double_dot(:)
real(dp)                       :: new_q(size(old_q))

if (present(old_q_double_dot)) then
   new_q(:) = old_q(:) + dT*old_q_dot(:) + dt**2*old_q_double_dot(:)/2._dp
else 
   new_q(:) = old_q(:) + dT*old_q_dot(:)
end if

end function

end module solver_utils
