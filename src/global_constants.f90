! ***********************************************************************************
! Module that contains all the constants whose values do not change during execution
!
! Constants can be floats such as PI
! Constant can be integers such as the number of bodies, number of joints
! ***********************************************************************************
module global_constants

  implicit none

  ! ****************************
  !  Define integer constants
  ! ****************************

  integer, parameter      :: SP = kind(0.0)    ! single precision
  integer(sp), parameter  :: DP = kind(0.0d0)  ! double precision

  integer(sp), parameter  :: NUM_SPAT_DIM = 3

  integer(sp), parameter  :: PRINT_LEVEL = 1 
  integer(sp), parameter  :: NUM_BODIES = 1
  integer(sp), parameter  :: NUM_JOINTS = 1


  integer(sp), parameter  :: MAX_TIME_STEPS = 100 ! for storing states

  integer(sp), parameter  :: NUM_DYNAM_EQN = 4 ! 2 kinematics and 2 dynamics in vector form per body
  integer(sp), parameter  :: NUM_ELAST_EQN = 0 ! number of elastic equatins in vector form per body
  integer(sp), parameter  :: NUM_JOINT_EQN = 0 ! number of joint equations in vector form per joint

  ! number of governing equations in vector form
  integer(sp), parameter  :: NUM_GOV_EQN = NUM_DYNAM_EQN + NUM_ELAST_EQN + NUM_JOINT_EQN ! per body

  ! number of state variables  
  integer(sp), parameter  :: NUM_STATES_PER_BODY = NUM_GOV_EQN*NUM_SPAT_DIM
  integer(sp), parameter  :: NUM_STATES = NUM_STATES_PER_BODY*NUM_BODIES
  
  ! ****************************
  ! Define real  constants
  ! ****************************

  real(dp), parameter     :: PI = 22.0_dp/7.0_dp
  real(dp), parameter     :: DEG_IN_RAD = 180.0_dp/PI
  real(dp), parameter     :: RAD_IN_DEG = PI/180.0_dp
  real(dp), parameter     :: ZERO = 0.0_dp

end module global_constants
