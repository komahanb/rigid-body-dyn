!*******************************************************************
! Module that helps share global variables across compilation units
!*******************************************************************
module global_variables
  use system_class
  use global_constants, only: MAX_TIME_STEPS, TOT_NDOF, ZERO
  use dispmodule, only: disp
  use tictoc, only: timer, timer_start, timer_stop
  use types, only: sp, dp, vector, matrix

  implicit none

  save

  ! stores the system in this object
  class(system), allocatable :: dynsys

  ! count the number of function and jacobian calls made
  integer(sp)  :: fcnt=0, fgcnt=0

  integer(sp)  :: MAX_NEWTON_ITER = 50

  ! acceleration due to gravity (inertial frame)
  type(vector), parameter :: GRAV = vector((/ 0.0_dp, -9.8_dp, 0.0_dp/))

  ! unit vector
  type(vector), parameter :: unitV = vector((/ ZERO, ZERO, ZERO /))

  ! zero vector
  type(vector), parameter :: zeroV = vector((/ ZERO, ZERO, ZERO /)) 

  logical  :: initialized = .false.

  ! integration time step and other times (these are just defaults)
  real(dp) :: dT         = 0.01_dp
  real(dp) :: time       = 0.0_dp
  real(dp) :: start_time = 0.0_dp
  real(dp) :: end_time   = 1.0_dp

  ! multiplicative factor in assembling the jacobian
  real(dp) :: aa != 1.0_dp/dT     ! default value
  real(dp) :: bb != 1.0_dp/dT**2  ! default value

  ! time history of states and their time derivative
  real(dp), dimension(TOT_NDOF, MAX_TIME_STEPS) :: q_save     = 0.0_dp
  real(dp), dimension(TOT_NDOF, MAX_TIME_STEPS) :: q_dot_save = 0.0_dp
  ! time history of residual
  real(dp), dimension(TOT_NDOF, MAX_TIME_STEPS) :: res_save   = 0.0_dp
  ! time history of residual update
  real(dp), dimension(TOT_NDOF, MAX_TIME_STEPS) :: dq_save    = 0.0_dp
  ! time history of jacobian
  real(dp), dimension(TOT_NDOF, TOT_NDOF, MAX_TIME_STEPS) ::  &
       &jac_save = 0.0_dp 

  !-------------------------------------------------------------------!
  ! Timing objects
  !-------------------------------------------------------------------!

  type(timer) :: time_jac    ! timer for jacobian assembly
  type(timer) :: time_res    ! timer for residual assemble
  type(timer) :: time_int    ! time for time integration with 
  type(timer) :: time_sol    ! time for each steady soluion=
  type(timer) :: time_tot    ! timer for whole calculation

  !-------------------------------------------------------------------!
  ! Error code
  !-------------------------------------------------------------------!
  integer :: ierr

  !-------------------------------------------------------------------!
  ! MPI parameter
  !-------------------------------------------------------------------!
  logical :: master = .false. ! am i master
  integer :: rank             ! rank of processor
  integer :: nproc            ! number of processors

  !-------------------------------------------------------------------!  
  ! solver type and other user settings
  !-------------------------------------------------------------------!

  character(len=25) :: solver_type

  logical           :: gradient    = .false.
  logical           :: fdgradient  = .false.
  logical           :: csgradient  = .false.
  logical           :: adjgradient = .false.

  
  !file handle
  integer(sp )                       :: filenum    

  ! output filename
  character                          :: filename   


contains

  !*******************************************************************!
  ! Performs initialization tasks such as loading input files and 
  ! sanity checks on the user supplied settings in global_constants
  !*******************************************************************!

  subroutine initialize()

    use mpi

    integer            :: i         ! loop counter
    integer            :: argc      ! number of command line arguments
    character(len=25)  :: argv       ! a command line argument

    !-----------------------------------------------------------------!
    ! initialize mpi
    !-----------------------------------------------------------------!

    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world,rank,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)

    if (ierr.ne.0) call finalize()

    ! find master
    if (rank .eq. 0) then       
       master = .true.
       call disp('>> Performing initialization...')
       call disp('>> Number of processors :', nproc)
    end if

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
    call mpi_finalize(ierr)

    if(master)  call disp (">> End of execution...")

  end subroutine finalize

  !*******************************************************************!
  ! Prints a total summary of the elapsed time for different opertns
  !************* *****************************************************!
  subroutine print_timer_summary()

    use dispmodule, only: disp

    !    call disp(">>------------------------------------------<<")
    call disp(">>-------- Timer Summary -------------------<<")
    !    call disp(">>------------------------------------------<<")

    ! time taken for residual assembly
    !    if (time_res%running) then
    call disp(" > Residual Assembly : ", (/ time_res%elapsed,&
         & time_res%elapsed/time_tot%elapsed /), orient="ROW")
    !    end if

    ! time taken for jacobian assembly
    !    if (time_jac%running) then
    call disp(" > Jacobian Assembly : ",  (/ time_jac%elapsed,&
         & time_jac%elapsed/time_tot%elapsed /), orient="ROW")
    !    end if

    ! time taken for Newton solution to steady state
    !    if (time_sol%running) then
    call disp(" > Newton Solution   : ",  (/ time_sol%elapsed,&
         & time_sol%elapsed/time_tot%elapsed /), orient="ROW")
    !    end if

    ! time taken for time marching
    !    if (time_int%running) then
    call disp(" > Time Integration  : ",  (/ time_int%elapsed,&
         & time_int%elapsed/time_tot%elapsed /), orient="ROW")
    !    end if

    ! time taken for total time
    !    if (time_tot%running) then
    call disp(" > Total Time        : ", (/ time_tot%elapsed,&
         & time_tot%elapsed/time_tot%elapsed /), orient="ROW")
    !    end if

    call disp(">>------------------------------------------<<")

  end subroutine print_timer_summary

end module global_variables
