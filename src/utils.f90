module utils
use constants

implicit none

contains

!--------------------------------------------------------
! returns a skew-symmetric matrix for doing cross product
! a  = [a1, a2, a3];
! \tilde{a}  = [0, -a_z, a_y; a_z, 0, -a_x; -a_y, a_x 0]
!--------------------------------------------------------
function skew(a) result(tilde_a)
  real(dp), intent(in)           :: a(3)
  real(dp)                       :: tilde_a(3,3)
  tilde_a=reshape((/ 0.0_dp, a(3), -a(2), -a(3), 0.0_dp, a(1),  a(2), -a(1), 0.0_dp /), (/3, 3/))
end function skew

!------------------------------------------
! returns cross product of vectors a and b
!------------------------------------------
function cross(a ,b) 
  real(dp), intent(in) :: a(3), b(3)
  real(dp)             :: cross(3)
  cross =matmul(skew(a), b) 
end function cross

end module utils
