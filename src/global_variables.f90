!*******************************************************************
! Module that helps share global variables across compilation units
! 
! 
!*******************************************************************
module global_variables

  use global_constants
  use types, only: sp, dp, vector, matrix
  use dispmodule, only: disp

  implicit none

  type(vector), parameter  :: GRAV = vector((/ ZERO, -1.0d0, ZERO /)) ! acceleration due to gravity
  type(vector), parameter :: unitV = vector((/ ZERO, ZERO, ZERO /)) ! unit vector
  type(vector), parameter :: zeroV = vector((/ ZERO, ZERO, ZERO /)) ! zero vector

  logical  :: elastic = .false.
  logical  :: initialized = .false.

  real(dp), dimension(NUM_STATES, MAX_TIME_STEPS) :: q_save=0.0d0, qdot_save=0.0d0 ! time history of states and their time derivative

  real(dp) :: aa = 1.0d0, bb =1.0d0 ! multiplicative factor in assembling te jacobian

  real(dp) :: dT ! integration time step

contains

  subroutine init()

    call disp('>> Performing initialization...')

    ! load input files may be?

    initialized = .true.

    call disp('>> Initialization complete...')

  end subroutine init

end module global_variables
