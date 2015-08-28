module solver_utils
use settings
implicit none

private
public get_updated_q, get_updated_q_dot, get_updated_q_double_dot,&
     & get_extrapolated_q

contains

!----------------------------------------------------
! updates the q vector value with the computed update
!----------------------------------------------------
function get_updated_q(old_q, del_q) result(new_q)
real(dp), intent(in)           :: del_q(:)
real(dp), intent(in)           :: old_q(:)
real(dp)                       :: new_q(size(old_q))

new_q(:) = old_q(:) + del_q(:)

end function

!--------------------------------------------------------
! updates the q_dot vector value with the computed update
!--------------------------------------------------------
function get_updated_q_dot(del_T, old_q_dot, del_q) result(new_q_dot)
real(dp), intent(in)           :: del_q(:)
real(dp), intent(in)           :: old_q_dot(:)
real(dp), intent(in)           :: del_T

real(dp)                       :: alpha0
real(dp)                       :: new_q_dot(size(old_q_dot))

new_q_dot(:) = old_q_dot(:) + alpha0*del_q(:)/del_t

end function
!---------------------------------------------------------------
! updates the q_double_dot vector value with the computed update
!---------------------------------------------------------------
function get_updated_q_double_dot(del_T, old_q_double_dot, del_q) result(new_q_double_dot)
real(dp), intent(in)           :: del_q(:)
real(dp), intent(in)           :: old_q_double_dot(:)
real(dp), intent(in)           :: del_T

real(dp)                       :: beta0
real(dp)                       :: new_q_double_dot(size(old_q_double_dot))

new_q_double_dot(:) = old_q_double_dot(:) + beta0*del_q(:)/del_t**2

end function

!--------------------------------------------------------------------------
! returns the extrapolated value of q based on first and second derivatives
!-------------------------------------------------------------------------
function get_extrapolated_q(del_T, old_q, old_q_dot, old_q_double_dot) result(new_q)
real(dp), intent(in)           :: del_T
real(dp), intent(in)           :: old_q(:), old_q_dot(:)
real(dp), intent(in), optional :: old_q_double_dot(:)
real(dp)                       :: new_q(size(old_q))

if (present(old_q_double_dot)) then
   new_q(:) = old_q(:) + del_T*old_q_dot(:) + del_T**2*old_q_double_dot(:)/2._dp
else 
   new_q(:) = old_q(:) + del_T*old_q_dot(:)
end if

end function

end module solver_utils
