module utils
use constants

implicit none

contains

!--------------------------------------------------------
! returns a skew-symmetric matrix for doing cross product
! a  = [a1, a2, a3];
! \tilde{a}  = [0, -a_z, a_y; a_z, 0, -a_x; -a_y, a_x 0]
!--------------------------------------------------------
function get_skew_sym_mat(a) result(tilde_a)
real(dp), intent(in)           :: a(3)
real(dp)                       :: tilde_a(3,3)

if (size(a) .ne. 3) then
   print*, "size of (a):",size(a)
   stop"wrong dim in the cross product"
end if

tilde_a=reshape((/ 0.0_dp, a(3), -a(2), -a(3), 0.0_dp, a(1),  a(2), -a(1), 0.0_dp /), (/3, 3/))

end function get_skew_sym_mat

end module utils
