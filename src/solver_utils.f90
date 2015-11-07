module solver_utils

use global_constants
use global_variables
!use utils

implicit none

private                                                               ! makes all functions default by private

public get_updated_q, get_updated_q_dot, get_updated_q_double_dot,&   ! expose only the functions that are needed outside the module
     & get_extrapolated_q, get_bdf_coeffs , get_approximated_q_dot,&
     & get_approximated_q_double_dot

contains

!----------------------------------------------------
! updates the q vector value with the computed update
!----------------------------------------------------
function get_updated_q(old_q, del_q) result(new_q)

real(dp), intent(in)           :: del_q(NUM_STATES)
real(dp), intent(in)           :: old_q(NUM_STATES)
real(dp)                       :: new_q(NUM_STATES)

new_q(:) = old_q(:) + del_q(:)

end function get_updated_q

!--------------------------------------------------------
! updates the q_dot vector value with the computed update
!--------------------------------------------------------
function get_updated_q_dot(old_q_dot, del_q) result(new_q_dot)

real(dp), intent(in)           :: del_q(NUM_STATES)
real(dp), intent(in)           :: old_q_dot(NUM_STATES)
real(dp)                       :: new_q_dot(size(old_q_dot))

new_q_dot(:) = old_q_dot(:) + aa*del_q(:)/dT

end function get_updated_q_dot

!---------------------------------------------------------------
! updates the q_double_dot vector value with the computed update
!---------------------------------------------------------------
function get_updated_q_double_dot(old_q_double_dot, del_q) result(new_q_double_dot)
real(dp), intent(in)           :: del_q(NUM_STATES)
real(dp), intent(in)           :: old_q_double_dot(NUM_STATES)
real(dp)                       :: new_q_double_dot(NUM_STATES)

new_q_double_dot(:) = old_q_double_dot(:) + bb*del_q(:)/dT**2

end function get_updated_q_double_dot


!--------------------------------------------------------------------------
! returns the extrapolated value of q based on first and second derivatives
!-------------------------------------------------------------------------
function get_extrapolated_q(old_q, old_q_dot, old_q_double_dot) result(new_q)
real(dp), intent(in)           :: old_q(NUM_STATES), old_q_dot(NUM_STATES)
real(dp), intent(in), optional :: old_q_double_dot(NUM_STATES)
real(dp)                       :: new_q(NUM_STATES)

if (present(old_q_double_dot)) then
   new_q(:) = old_q(:) + dT*old_q_dot(:) + dT**2*old_q_double_dot(:)/2._dp
else 
   new_q(:) = old_q(:) + dT*old_q_dot(:)
end if

end function get_extrapolated_q

!------------------------------------------------------------------
! returns the bdf coeffs (unscaled with respect to the step size h)
!------------------------------------------------------------------
function get_bdf_coeffs(d, m) result(c)
integer(sp), intent(in)    :: d             ! d-th derivative e.g. first or second derivative
integer(sp), intent(in)    :: m             ! m-th order accurate
integer(sp)                :: n             ! number of points needed for the required degree and accuracy
real(dp)                   :: c(d+m)        ! vector of coefficients used to approximate derivate
real(dp)                   :: x(d+m)        ! vector of evaluation points (not needed here)
real(dp), parameter        :: h = 1._dp     ! if we set h=dT we will get the scaled coeffs
integer(sp)                :: i

n = m +d                                    ! number of points needed for the reqd accuracy and degree
call differ_backward ( h, d, m, c, x )      ! calling a library function
!$!c = reverse_real(c)                         ! store in backward order for convenience
!!$x = reverse_real(x)                         ! store in backward order for convenience

end function get_bdf_coeffs

!# tested OK, just need to check the sign of the derivative
!# how to deal with initial steps?
!--------------------------------------------------------
! returns the approximated first derivative
! use 'q' to produce an m-th order approximation to q_dot
!--------------------------------------------------------
function get_approximated_q_dot(q, m) result(q_dot)
!************************************************
! q = [t= 0                , (q1, q2, q_{NUM_STATES}), 
!      t= 0+dT          , (q1, q2, q_{NUM_STATES}),
!         .                ,          .         ,
!         .                ,          .         ,
!      t= 0+k*dT        , (q1, q2, q_{NUM_STATES})]
!************************************************
integer(sp), parameter   :: degree = 1           ! since we are approximating first derivative
integer(sp)              :: cnt                  ! cnt = m + degree -1
integer(sp), intent(in)  :: m                    ! m = order of accuracy of sought derivative
real(dp), intent(in)     :: q(0:degree+m-1,NUM_STATES)! matrix whose structure is drawn above
real(dp)                 :: q_dot(NUM_STATES)         ! output first derivative vector
real(dp)                 :: alpha(0:m+degree-1)  ! should always be based on the reqd. accuracy (not on the total available data like alpha(0 to k))
integer(sp)              :: i

! getting the BDF coefficients
alpha = get_bdf_coeffs(degree,m)

cnt = m + degree -1
if (size(alpha).ne. cnt+1) stop"Wrong operation predicted. stopping" ! use something diff?

q_dot(:)=0._dp                                   ! initialize
do i = 0, cnt                                    ! loop (sum) across the data points (0 to m in paper) 
   q_dot(:) =  q_dot(:) + alpha(i)*q(i,:)        ! find the cumulative sum
end do
q_dot(:) = q_dot(:)/dT

end function get_approximated_q_dot

!# tested OK, just need to check the sign of the derivative
!# how to deal with initial steps?
!--------------------------------------------------------
! returns the approximated second derivative
! use 'q' to produce an m-th order approximation to q_double_dot
!--------------------------------------------------------
function get_approximated_q_double_dot(q, m) result(q_double_dot)
!------------------------------------------------------
! q = [t= 0                , (q1, q2, q_{NUM_STATES}), 
!      t= 0+dT          , (q1, q2, q_{NUM_STATES}),
!         .                ,          .         ,
!         .                ,          .         ,
!      t= 0+k*dT        , (q1, q2, q_{NUM_STATES})]
!-----------------------------------------------------
integer(sp), parameter   :: degree = 2           ! since we are approximating second derivative
integer(sp)              :: cnt                  ! cnt = m + degree -1
integer(sp), intent(in)  :: m                    ! m = order of accuracy of sought derivative
real(dp), intent(in)     :: q(0:degree+m-1,NUM_STATES)! matrix whose structure is drawn above
real(dp)                 :: q_double_dot(NUM_STATES)  ! output second derivative vector
real(dp)                 :: beta(0:m+degree-1)   ! should always be based on the reqd. accuracy (not on the total available data like beta(0 to k))
integer(sp)              :: i

! getting the BDF coefficients
beta = get_bdf_coeffs(degree,m)

cnt = m + degree -1
if (size(beta).ne. cnt+1) stop"Wrong operation predicted. stopping" ! use something diff?

q_double_dot(:)=0._dp                                         ! initialize
do i = 0, cnt                                                 ! loop (sum) across the data points (0 to m in paper) 
   q_double_dot(:) =  q_double_dot(:) + beta(i)*q(i,:)        ! find the cumulative sum
end do
q_double_dot(:) = q_double_dot(:)/dT**2

end function get_approximated_q_double_dot


end module solver_utils
