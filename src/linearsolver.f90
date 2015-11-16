module linearsolver 

!-module options

  implicit none

contains

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

    ! create temp working vectors for GMRES (see dslugm.f for details)
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

end module linearsolver 
