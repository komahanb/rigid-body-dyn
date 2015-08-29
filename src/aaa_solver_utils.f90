module solver_utils
use settings
implicit none

private
public get_updated_q, get_updated_q_dot, get_updated_q_double_dot,&
     & get_extrapolated_q, get_bdf_coeffs

contains

!----------------------------------------------------
! updates the q vector value with the computed update
!----------------------------------------------------
function get_updated_q(old_q, del_q) result(new_q)
real(dp), intent(in)           :: del_q(:)
real(dp), intent(in)           :: old_q(:)
real(dp)                       :: new_q(size(old_q))

new_q(:) = old_q(:) + del_q(:)

end function get_updated_q

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

end function get_updated_q_dot
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

end function get_updated_q_double_dot

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

end function get_extrapolated_q

!------------------------------------------------------------------
! returns the bdf coeffs (unscaled with respect to the step size h)
!------------------------------------------------------------------
function get_bdf_coeffs(d, m) result(c)

integer(sp), intent(in)    :: d           ! d-th derivative e.g. first or second derivative
integer(sp), intent(in)    :: m           ! m-th order accurate
!integer(sp)     :: n = d + m   ! number of points needed for the required degree and accuracy
real(dp)                   :: c(d+m)        ! vector of coefficients used to approximate derivate
real(dp)                   :: x(d+m)        ! vector of evaluation points (not needed here)
real(dp), parameter        :: h = 1._dp   ! if we set h=delta_T we will get the scaled coeffs
integer(sp)                :: i

call differ_backward ( h, d, m, c, x )    ! calling a library function
write (*, '(a,i2,a,i2)' )  '  Backward difference coefficients, d = ', d, ', m = ', m
write (*, *) "index  ",  "coeff"
do i = 1, d + m
   write (*, '(a,i4,f13.2)') 'k', int(x(i)) , c(i)
end do

end function get_bdf_coeffs


end module solver_utils
