program main

  ! program references

  !  use global_constants, only: dp
  !  use global_variables, only: nproc, rank
  !  use jacobian_matrix_class, only: jacobian_matrix
  !  use matrix_class, only: matrix
  !  use tictoc, only: timer_start, timer_stop, timer
  use global_variables, only: print_timer_summary, master, nproc, rank,&
       & initialized, solver_type, time_tot
  use tictoc, only: timer_start, timer_stop
  use dispmodule, only: disp

  implicit none
#include <finclude/petsc.h90>
!#include <finclude/slepcsys.h>
!#include <finclude/slepceps.h>

!#include <finclude/petsc.h90>
!#include <slepc/finclude/slepcsys.h>
!#include <slepc/finclude/slepceps.h>
!#include <finclude/petscdef.h>

  Mat         :: jac         ! jacobian matrix
  ! Mat         :: jac_prec    ! preconditioned jacobian matrix
  Vec         :: resvec      ! residual vector
  Vec         :: xvec        ! results
  KSP         :: ksp         ! linear solver context
  PC          :: pc          ! preconditioner
  SNES        :: snes        ! nonlinear solver context
  integer     :: ierr        ! error flag

  !-------------------------------------------------------------------!
  ! initialize
  !-------------------------------------------------------------------!

  call initialize()

!  call PetscOptionsSetValue("-snes_mf_operator","TRUE",ierr)

  !-------------------------------------------------------------------!
  ! Execute multibody dynamics
  !-------------------------------------------------------------------!

  ! call execute_dynamics()

  !-------------------------------------------------------------------!
  ! finalize run
  !-------------------------------------------------------------------!

  call finalize()


contains


 subroutine initialize()

!    use mpi

    integer            :: i         ! loop counter
    integer            :: argc      ! number of command line arguments
    character(len=25)  :: argv ,cc      ! a command line argument

    !-----------------------------------------------------------------!
    ! initialize mpi
    !-----------------------------------------------------------------!

    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

    call MPI_COMM_RANK(PETSC_COMM_WORLD,rank,ierr)
    call MPI_COMM_SIZE(PETSC_COMM_WORLD,nproc,ierr)

!   call mpi_init(ierr)
!   call mpi_comm_rank(mpi_comm_world,rank,ierr)
!   call mpi_comm_size(mpi_comm_world,nproc,ierr)
    
    if (ierr.ne.0) call finalize()

    ! find master
    if (rank .eq. 0) then       
       master = .true.
       call disp('>> Performing initialization...')
       call disp('>> Number of processors :', nproc)
    end if
    print*, "I am processor", rank , " of",nproc

    !-----------------------------------------------------------------!
    ! Start the total timer
    !-----------------------------------------------------------------!
    call timer_start(time_tot)
!
    !-----------------------------------------------------------------!
    ! process command line arguments
    !-----------------------------------------------------------------!

    ! get number of command line arguments
    argc = command_argument_count()

    ! begin loop around command line args
    do i = 1, argc

       ! get that argument
       call get_command_argument(i,argv)

       ! begin case structure
       select case(trim(argv))

          ! solver
       case('--solver_type')

          ! get next argument
          call get_command_argument(i+1,argv)

          ! set global var
          solver_type = trim(argv)

          ! do nothing here
       case default

       end select

    end do

    initialized = .true.

    if (master) call disp('>> Initialization complete...')

  end subroutine initialize


  !*******************************************************************!
  ! finalize
  !*******************************************************************!

  subroutine finalize()

    !-----------------------------------------------------------------!
    ! Stop the timer
    !-----------------------------------------------------------------!
    if (master) then
       call timer_stop(time_tot)
       call print_timer_summary()
    end if

    !-----------------------------------------------------------------!
    ! finalize MPI
    !-----------------------------------------------------------------!
    ! call mpi_finalize(ierr)
    
    call PetscFinalize(PETSC_NULL_CHARACTER,ierr)
    
    if(master)  call disp (">> End of execution...")

  end subroutine finalize


end program main
