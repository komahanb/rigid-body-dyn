!==============================================================================
! Module that solves the linear system AX=b using different third-party solvers
!
! You may want to write interfaces to solvers here
!==============================================================================
module linear_system

  use global_constants

  implicit none

  !make all solvers by default private
  !  private 

  !  public

contains

  ! ===========================================================================
  ! Returns the inverse of a matrix calculated by finding the LU decomposition
  ! ===========================================================================
  function inv(A) result(Ainv)

    real(dp), dimension(:,:), intent(in)     :: A      ! A matrix
    real(dp), dimension(size(A,1),size(A,2)) :: Ainv
    real(dp), dimension(size(A,1))           :: work   ! work array for LAPACK
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

  ! ==========================================================================
  ! Returns the inverse of a matrix calculated by finding the SV decomposition
  ! ==========================================================================
  function inv2(A,b,n) result(x)
    implicit none
    integer n
    real(dp) :: x(n),b(n)
    real(dp) :: A(n,n)
    real(dp) :: err
    integer  :: i, info, lda, ldb, nrhs
    integer  :: ipiv(n)

    external DGESV

    nrhs = 1 ! number of right hand sides in b
    lda = n  ! leading dimension of a
    ldb = n  ! leading dimension of b

    call dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)

    if( info.gt.0 ) then
       write(*,*)'The diagonal element of the triangular factor of a,'
       write(*,*)'u(',info,',',info,') is zero, so that'
       write(*,*)'a is singular; the solution could not be computed.'
       stop
    end if

    ! Note: the solution is returned in b  and a has been changed.
    x=b

  end function inv2

end module linear_system
