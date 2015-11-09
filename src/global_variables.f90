!*******************************************************************
! Module that helps share global variables across compilation units
!*******************************************************************
module global_variables

  use global_constants
  use types, only: sp, dp, vector, matrix
  use dispmodule, only: disp

  implicit none

  integer(sp)  :: MAX_NEWTON_ITER = 50

  ! acceleration due to gravity (inertial frame)
  type(vector), parameter :: GRAV = vector((/ ZERO, -1.0_dp, ZERO /))

  ! unit vector
  type(vector), parameter :: unitV = vector((/ ZERO, ZERO, ZERO /))

  ! zero vector
  type(vector), parameter :: zeroV = vector((/ ZERO, ZERO, ZERO /)) 

  logical  :: initialized = .false.

  ! integration time step and other times (these are just defaults)
  real(dp) :: dT         = 0.01_dp
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

contains
  
  !*******************************************************************!
  ! Performs initialization tasks such as loading input files and 
  ! sanity checks on the user supplied settings in global_constants
  !*******************************************************************!
  subroutine init()

    call disp('>> Performing initialization...')

    ! load input files may be and make appropriate settings

    ! implement logic to check if the user entered settings are good

    initialized = .true.

    call disp('>> Initialization complete...')

  end subroutine init

end module global_variables
