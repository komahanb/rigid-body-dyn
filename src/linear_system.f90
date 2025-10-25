!==============================================================================
! Module that solves the linear system AX=b using different third-party solvers
!
! You may want to write interfaces to solvers here
!==============================================================================
module linear_system

  use global_constants

  implicit none

  !make all solvers by default private
! private 

!  public :: solve

  interface direct_solve
     module procedure inv, inv2
  end interface direct_solve


  !interface iter_solve
  !   module procedure gmres_solve
  !end interface iter_solve

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
    integer  :: info, lda, ldb, nrhs
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

  !===================================================================!
  ! GMRES wrapper for DSLAP
  !===================================================================!
  subroutine gmres(n, nnz, row, col, val, x, b, tol, err, maxiter, iter)

    !---arguments

    integer, intent(in)    :: n        ! order of matrix
    integer, intent(in)    :: nnz      ! number of nonzeros in sparse matrix
    integer, intent(in)    :: row(nnz) ! vector of rows
    integer, intent(in)    :: col(nnz) ! vector of columns
    integer, intent(in)    :: maxiter  ! max iterations
    integer, intent(out)   :: iter     ! actual iterations

    real(8), intent(in)    :: val(nnz) ! vector of matrix values
    real(8), intent(inout) :: x(n)     ! unknown vector
    real(8), intent(inout) :: b(n)     ! right hand side vector
    real(8), intent(out)   :: err      ! actual error on output
    real(8), intent(in)    :: tol      ! linear inner tolerance

    !---local variables

    integer :: ierr                  ! error code
    integer :: lenw                  ! length of rwork
    integer :: leniw                 ! length of iwork
    integer :: nsave = 20            ! vectors to orthogonalize with
    integer, allocatable :: iwork(:) ! work vector for integer
    real(8), allocatable :: rwork(:) ! working vector for real

    !---begin execution

    ! create temp working vectors for GMRES 
    ! (see dslugm.f for details)
    lenw =  1 + n*(nsave+7) + nsave*(nsave+3) + nnz 
    leniw = 32 + 4*n + nnz

    allocate(rwork(lenw))
    allocate(iwork(leniw))

    ! solve matrix with GMRES
    !call DSLUGM(n, b, x, nnz, row, col, val, 0, nsave, 0, tol, maxiter, iter, &
    !     err, ierr, 0, rwork, lenw, iwork, leniw)

    ! return allocated memory
    deallocate(rwork)
    deallocate(iwork)

  end subroutine gmres

!!$  ! Call the GMRES solver to calculate the state update
!!$  subroutine gmres_solve( )
!!$
!!$    !---local variables
!!$
!!$    integer :: i ! iteration index
!!$    integer :: maxiter
!!$    integer :: iter
!!$    real(8) :: errork
!!$    real(8) :: errorf
!!$    real(8) :: err
!!$    real(8) :: knew ! new eigenvalue
!!$    real(8) :: kold ! old eigenvalue
!!$    real(8), allocatable :: qnew(:) ! new state
!!$    real(8), allocatable :: qold(:) ! old state

    !    real(8), allocatable :: sold(:) ! old sourcie
    !    real(8), allocatable :: snew(:) ! new source
    
    ! convert sparse matrix notation for linear solver
!!$    call DS2Y(loss_matrix % n, loss_matrix % nz, &
!!$         &loss_matrix % row, &
!!$         loss_matrix % col, loss_matrix % val, 0)
!!$
!!$    call DS2Y(prod_matrix % n, prod_matrix % nz, &
!!$         &prod_matrix % row, &
!!$         prod_matrix % col, prod_matrix % val, 0)
    
    ! allocate temporary variables
!!$    allocate(qnew(TOT_NDOF))
!!$    allocate(qold(TOT_NDOF))
!!$
!!$    maxiter = 100
!!$
!!$    ! solve for new flux estimate using gmres linear solver
!!$    call gmres(TOT_NDOF,  sds, loss_matrix % row, &
!!$         loss_matrix % col, loss_matrix % val, qnew, sold, &
!!$         itol, err, maxiter, iter)
!!$  
!!$    ! calculate new fission source
!!$    call DSMV(prod_matrix % n, phinew, snew, prod_matrix % nz, &
!!$         prod_matrix % row, prod_matrix % col, &
!!$         prod_matrix % val, 0)
!!$
!!$    deallocate(qnew, qold)
!!$
!!$  end subroutine solve_system


end module linear_system
