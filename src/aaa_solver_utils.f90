module solver_utils
use constants
use common_utils
implicit none

private                                                               ! makes all functions default by private

public get_updated_q, get_updated_q_dot, get_updated_q_double_dot,&   ! expose only the functions that are needed outside the module
     & get_extrapolated_q, get_bdf_coeffs !, get_approximated_q_dot

contains

!# yet to be tested
!----------------------------------------------------
! updates the q vector value with the computed update
!----------------------------------------------------
function get_updated_q(old_q, del_q) result(new_q)
real(dp), intent(in)           :: del_q(dim_q)
real(dp), intent(in)           :: old_q(dim_q)
real(dp)                       :: new_q(dim_q)

new_q(:) = old_q(:) + del_q(:)

end function get_updated_q

!# yet to be tested
!--------------------------------------------------------
! updates the q_dot vector value with the computed update
!--------------------------------------------------------
function get_updated_q_dot(old_q_dot, del_q) result(new_q_dot)
real(dp), intent(in)           :: del_q(dim_q)
real(dp), intent(in)           :: old_q_dot(dim_q)

real(dp), parameter            :: alpha0 = 1._dp
real(dp)                       :: new_q_dot(size(old_q_dot))

new_q_dot(:) = old_q_dot(:) + alpha0*del_q(:)/del_t

end function get_updated_q_dot

!# yet to be tested
!---------------------------------------------------------------
! updates the q_double_dot vector value with the computed update
!---------------------------------------------------------------
function get_updated_q_double_dot(old_q_double_dot, del_q) result(new_q_double_dot)
real(dp), intent(in)           :: del_q(dim_q)
real(dp), intent(in)           :: old_q_double_dot(dim_q)
real(dp), parameter            :: beta0 = 1._dp
real(dp)                       :: new_q_double_dot(dim_q)

new_q_double_dot(:) = old_q_double_dot(:) + beta0*del_q(:)/del_t**2

end function get_updated_q_double_dot

!# tested OK
!--------------------------------------------------------------------------
! returns the extrapolated value of q based on first and second derivatives
!-------------------------------------------------------------------------
function get_extrapolated_q(old_q, old_q_dot, old_q_double_dot) result(new_q)
real(dp), intent(in)           :: old_q(dim_q), old_q_dot(dim_q)
real(dp), intent(in), optional :: old_q_double_dot(dim_q)
real(dp)                       :: new_q(dim_q)

if (present(old_q_double_dot)) then
   new_q(:) = old_q(:) + del_t*old_q_dot(:) + del_t**2*old_q_double_dot(:)/2._dp
else 
   new_q(:) = old_q(:) + del_t*old_q_dot(:)
end if

end function get_extrapolated_q

!# yet to be tested
!------------------------------------------------------------------
! returns the bdf coeffs (unscaled with respect to the step size h)
!------------------------------------------------------------------
function get_bdf_coeffs(d, m) result(c)
integer(sp), intent(in)    :: d             ! d-th derivative e.g. first or second derivative
integer(sp), intent(in)    :: m             ! m-th order accurate
integer(sp)                :: n             ! number of points needed for the required degree and accuracy
real(dp)                   :: c(d+m)        ! vector of coefficients used to approximate derivate
real(dp)                   :: x(d+m)        ! vector of evaluation points (not needed here)
real(dp), parameter        :: h = 1._dp     ! if we set h=del_t we will get the scaled coeffs
integer(sp)                :: i

n = m +d                                    ! number of points needed for the reqd accuracy and degree

call differ_backward ( h, d, m, c, x )      ! calling a library function

c = reverse_real(c)                         ! store in backward order for convenience
!!$x = reverse_real(x)                         ! store in backward order for convenience
!!$write (*, '(a,i2,a,i2)' )  '  Backward difference coefficients, d = ', d, ', m = ', m
!!$write (*, *) "index  ",  "coeff"
!!$
!!$do i = 1, n
!!$   write (*, '(a,i4,f13.2)') 'k', int(x(i)) , c(i)
!!$end do

end function get_bdf_coeffs

!------------------------------------------
! returns the approximated first derivative
!------------------------------------------
!!$function get_approximated_q_dot(q, nvars, npts) result(q_dot)
!!$
!!$integer(sp), intent(in)  :: nvars, npts
!!$real(dp), intent(in)     :: q(nvars npts)
!!$integer(sp), parameter   :: degree = 1 
!!$real(dp)                 :: q_dot(nvars)
!!$real(dp), allocatable    :: alpha(npts)
!!$real(dp)                 :: total = 0._dp
!!$integer(sp)              :: i
!!$
!!$end function get_approximated_q_dot

end module solver_utils
