!=====================================================================!
! Module that contains helper functions used while solving time-
! integration problem. 
!---------------------------------------------------------------------!
! Contains functionst to:
! (a) compute the extrapolations of state based on previous values
! (b) update the state, its first and second time derivatives
! (c) get the coefficients for BDF
!=====================================================================!
module solver_utils

use global_constants, only: sp, dp, TOT_NDOF
use global_variables, only: dt, aa, bb
!use utils

implicit none

! makes all functions default by private
private                                                               

! expose only the functions that are needed outside the module
public get_updated_q, get_updated_q_dot, get_updated_q_double_dot,&   
     & get_approx_q, get_bdf_coeffs , get_approx_q_dot,&
     & get_approx_q_double_dot

contains

!----------------------------------------------------
! updates the q vector value with the computed update
!----------------------------------------------------
function get_updated_q(old_q, del_q) result(new_q)

real(dp), intent(in)           :: del_q(TOT_NDOF)
real(dp), intent(in)           :: old_q(TOT_NDOF)
real(dp)                       :: new_q(TOT_NDOF)

new_q(:) = old_q(:) + del_q(:)

end function get_updated_q

!--------------------------------------------------------
! updates the q_dot vector value with the computed update
!--------------------------------------------------------
function get_updated_q_dot(old_q_dot, del_q) result(new_q_dot)

real(dp), intent(in)           :: del_q(TOT_NDOF)
real(dp), intent(in)           :: old_q_dot(TOT_NDOF)
real(dp)                       :: new_q_dot(size(old_q_dot))

!aa=alpha0/dT (globally set)

new_q_dot(:) = old_q_dot(:) + aa*del_q(:)   

end function get_updated_q_dot

!---------------------------------------------------------------
! updates the q_double_dot vector value with the computed update
!---------------------------------------------------------------
function get_updated_q_double_dot(old_q_double_dot, del_q)&
     & result(new_q_double_dot)
real(dp), intent(in)           :: del_q(TOT_NDOF)
real(dp), intent(in)           :: old_q_double_dot(TOT_NDOF)
real(dp)                       :: new_q_double_dot(TOT_NDOF)

!bb=beta0/dT**2 (globally set)

new_q_double_dot(:) = old_q_double_dot(:) + bb*del_q(:)

end function get_updated_q_double_dot


!--------------------------------------------------------------------------
! returns the extrapolated value of q based on first and second derivatives
!-------------------------------------------------------------------------
function get_approx_q(old_q, old_q_dot, old_q_double_dot) &
     &result(new_q)
real(dp), intent(in)           :: old_q(TOT_NDOF), old_q_dot(TOT_NDOF)
real(dp), intent(in), optional :: old_q_double_dot(TOT_NDOF)
real(dp)                       :: new_q(TOT_NDOF)

if (present(old_q_double_dot)) then
   new_q(:) = old_q(:) + dT*old_q_dot(:) &
        &+ dT**2*old_q_double_dot(:)/2.0_dp
else 
   new_q(:) = old_q(:) + dT*old_q_dot(:)
end if

end function get_approx_q

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
function get_approx_q_dot(q, m) result(q_dot)
!************************************************
! q = [t= 0                , (q1, q2, q_{TOT_NDOF}), 
!      t= 0+dT          , (q1, q2, q_{TOT_NDOF}),
!         .                ,          .         ,
!         .                ,          .         ,
!      t= 0+k*dT        , (q1, q2, q_{TOT_NDOF})]
!************************************************
integer(sp), parameter   :: degree = 1           ! since we are approximating first derivative
integer(sp)              :: cnt                  ! cnt = m + degree -1
integer(sp), intent(in)  :: m                    ! m = order of accuracy of sought derivative
real(dp), intent(in)     :: q(0:degree+m-1,TOT_NDOF)! matrix whose structure is drawn above
real(dp)                 :: q_dot(TOT_NDOF)         ! output first derivative vector
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

end function get_approx_q_dot

!# tested OK, just need to check the sign of the derivative
!# how to deal with initial steps?
!--------------------------------------------------------
! returns the approximated second derivative
! use 'q' to produce an m-th order approximation to q_double_dot
!--------------------------------------------------------
function get_approx_q_double_dot(q, m) result(q_double_dot)
!------------------------------------------------------
! q = [t= 0                , (q1, q2, q_{TOT_NDOF}), 
!      t= 0+dT          , (q1, q2, q_{TOT_NDOF}),
!         .                ,          .         ,
!         .                ,          .         ,
!      t= 0+k*dT        , (q1, q2, q_{TOT_NDOF})]
!-----------------------------------------------------
integer(sp), parameter   :: degree = 2           ! since we are approximating second derivative
integer(sp)              :: cnt                  ! cnt = m + degree -1
integer(sp), intent(in)  :: m                    ! m = order of accuracy of sought derivative
real(dp), intent(in)     :: q(0:degree+m-1,TOT_NDOF)! matrix whose structure is drawn above
real(dp)                 :: q_double_dot(TOT_NDOF)  ! output second derivative vector
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

end function get_approx_q_double_dot


end module solver_utils
