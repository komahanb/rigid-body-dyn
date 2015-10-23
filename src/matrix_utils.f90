!=======================================================================
! Module that contains utitilies needed for matrix operations
!
! identity, zero and unit matrices
! rotation and translation matrices !?? derivatives and transposes
!
!=======================================================================
module matrix_utils
  use constants
  implicit none

contains

  !************************************************
  ! returns an identity matrix of size num_spat_dim
  !************************************************
  function idM()
    type(matrix) :: idM
    idM = matrix(eye(num_spat_dim))
  end function idM

  !************************************************
  ! returns an zero matrix of size num_spat_dim  
  !************************************************
  function zeroM()
    type(matrix) :: zeroM
    zeroM = matrix(zeros(num_spat_dim))
  end function zeroM

  !************************************************
  ! returns an unit matrix of size num_spat_dim
  !************************************************
  function unitM()
    type(matrix) :: unitM
    unitM =  matrix(ones(num_spat_dim))
  end function unitM
  !************************************************
  ! generates a nxn identity matrix
  !************************************************
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
  !************************************************
  ! generates a nxn zero matrix
  !************************************************
  function zeros(n)
    integer(sp) :: n
    real(dp)    :: zeros(n,n)
    zeros        = 0.0_dp
  end function zeros

  !************************************************
  ! generates a nxn unit matrix
  !************************************************
  function ones(n)
    integer(sp) :: n
    real(dp)    :: ones(n,n)
    ones        = 1.0_dp
  end function ones

  !************************************************
  ! rotates a vector from body frame inertial
  !************************************************
  function CBI(theta)
    real(dp)                :: theta(num_spat_dim)
    real(dp)                :: CBI(num_spat_dim, num_spat_dim)

    CBI(1,1) =  cos(theta(2))*cos(theta(3)) + sin(theta(1))*sin(theta(2))*sin(theta(3))
    CBI(2,1) = -cos(theta(2))*sin(theta(3)) + sin(theta(1))*sin(theta(2))*cos(theta(3))
    CBI(3,1) =  cos(theta(1))*sin(theta(2))

    CBI(1,2) =  cos(theta(1))*sin(theta(3))
    CBI(2,2) =  cos(theta(1))*cos(theta(3))
    CBI(3,2) = -sin(theta(1))

    CBI(1,3) = -sin(theta(2))*cos(theta(3)) + sin(theta(1))*cos(theta(2))*sin(theta(3))
    CBI(2,3) =  sin(theta(2))*sin(theta(3)) + sin(theta(1))*cos(theta(2))*cos(theta(3))
    CBI(3,3) =  cos(theta(1))*cos(theta(2))

  end function CBI

  !************************************************
  ! transformation matrix 
  ! angular velocity between intertial and body axis
  !************************************************
  function SIB(theta)

    real(dp), intent(in)    :: theta(num_spat_dim)
    real(dp)                :: SIB(num_spat_dim, num_spat_dim)

    SIB(1,1) = cos(theta(3))
    SIB(2,1) = -sin(theta(3))
    SIB(3,1) = 0.0_dp

    SIB(1,2) = cos(theta(1))*sin(theta(3))
    SIB(2,2) = cos(theta(1))*cos(theta(3))
    SIB(3,2) = -sin(theta(1)) 

    SIB(1,3) = 0.0_dp
    SIB(2,3) = 0.0_dp
    SIB(3,3) = 1.0_dp

  end function SIB


!!$ check 
  function SINV (theta)

    real(dp), intent(in)    :: theta(num_spat_dim)
    real(dp)                :: SINV(num_spat_dim, num_spat_dim)

    SINV =  inv(SIB(theta))

!!$
!!$    SINV(1,1) = -sin(theta(3)) *theta_dot !cos(theta(3))
!!$    SINV(2,1) = -cos(theta(3))*theta_dot
!!$    SINV(3,1) = 0.0_dp
!!$
!!$    SINV(1,2) = cos(theta(1))*sin(theta(3)) 
!!$    SINV(2,2) = cos(theta(1))*cos(theta(3))
!!$    SINV(3,2) = -sin(theta(1)) 
!!$   
!!$    SINV(1,3) = 0.0_dp
!!$    SINV(2,3) = 0.0_dp
!!$    SINV(3,3) = 0.0_dp
!!$

  end function SINV


  function SDOT(theta)

    real(dp), intent(in)    :: theta(num_spat_dim)
    real(dp)                :: SDOT(num_spat_dim, num_spat_dim)

    SDOT(1,1) = -sin(theta(3))
    SDOT(2,1) = -cos(theta(3))
    SDOT(3,1) = 0.0_dp

    SDOT(1,2) = -sin(theta(1))*sin(theta(3)) + cos(theta(1))*cos(theta(3))
    SDOT(2,2) = -sin(theta(1))*cos(theta(3)) - cos(theta(1))*sin(theta(3))
    SDOT(3,2) = -cos(theta(1)) 

    SDOT(1,3) = 0.0_dp
    SDOT(2,3) = 0.0_dp
    SDOT(3,3) = 0.0_dp

    stop
  end function SDOT

  ! Returns the inverse of a matrix calculated by finding the LU
  ! decomposition.  Depends on LAPACK.
  function inv(A) result(Ainv)

    real(dp), dimension(:,:), intent(in)     :: A
    real(dp), dimension(size(A,1),size(A,2)) :: Ainv

    real(dp), dimension(size(A,1))           :: work  ! work array for LAPACK
    integer(sp), dimension(size(A,1))        :: ipiv   ! pivot indices
    integer(sp)                              :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if
  end function inv

end module matrix_utils
