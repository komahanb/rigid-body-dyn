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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!$! identity matrix
!!$function eye1()
!!$  type(matrix)    :: eye1
!!$  integer         :: i,j
!!$  eye1%ij = 0.0_dp
!!$  do i = 1, num_spat_dim
!!$     do j = 1, num_spat_dim
!!$        if (j.eq.i) eye1%ij(j,i) = 1.0_dp
!!$     end do
!!$  end do
!!$end function eye1
!!$
!!$! generates a nxn zero matrix
!!$function zeros1()
!!$   type(matrix)    :: zeros1
!!$   zeros1%ij(:,:)        = 0.0_dp
!!$end function zeros1
!!$
!!$! generates a nxn unit matrix
!!$function ones1()
!!$  type(matrix) :: ones1
!!$  ones1%ij        = 1.0_dp
!!$end function ones1


end module matrix_utils
