module matrix_utils
use constants

contains

! generates a nxn identity matrix
function eye(n)
  integer(sp) :: n
  real(dp)    :: eye(n,n)
  integer     :: i, j
  eye = 0.0_dp
  do i = 1, n
     do j = 1, n
        if (j.eq.i) eye(j,i) = 1.0_dp
     end do
  end do
end function eye

! generates a nxn zero matrix
function zeros(n)
  integer(sp) :: n
  real(dp)    :: zeros(n,n)
  zeros        = 0.0_dp
end function zeros

! generates a nxn unit matrix
function ones(n)
  integer(sp) :: n
  real(dp)    :: ones(n,n)
  ones        = 1.0_dp
end function ones

end module matrix_utils
