!==============================================================================
! Module that solves the linear system AX=b using different third-party solvers
!
! You may want to write interfaces to solvers here
!==============================================================================
module linear_system

  use global_constants

  implicit none

  !make all solvers by default private
  private 

  public :: solve

  interface solve
     module procedure solve_system
  end interface solve

contains

  ! Call the GMRES solver to calculate the state update
  subroutine solve_system( )

    !---local variables

    integer :: i ! iteration index
    integer :: maxiter
    integer :: iter
    real(8) :: errork
    real(8) :: errorf
    real(8) :: err
    real(8) :: knew ! new eigenvalue
    real(8) :: kold ! old eigenvalue
    real(8), allocatable :: qnew(:) ! new state
    real(8), allocatable :: qold(:) ! old state

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
    allocate(qnew(TOT_NDOF))
    allocate(qold(TOT_NDOF))

    maxiter = 100
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

    deallocate(qnew, qold)

  end subroutine solve_system

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
    call DSLUGM(n, b, x, nnz, row, col, val, 0, nsave, 0, tol, maxiter, iter, &
         err, ierr, 0, rwork, lenw, iwork, leniw)

    ! return allocated memory
    deallocate(rwork)
    deallocate(iwork)

  end subroutine gmres


end module linear_system
