!=====================================================================!
! Module that contains all the constants whose values do not change 
! during execution
!---------------------------------------------------------------------!
! Example:
!
! Constants can be doubles such as PI
! Constant can be integers such as number of bodies, number of joints
!
! Usage: 
!
! TUNABLE params can be changed as per the needs
! NON-TUNABLE params can not be changed without appropriate code change
!=====================================================================!
module global_constants

  implicit none

  ! make everything private
  private
 
  ! expose double params
  public PI, DEG_PER_RAD, RAD_PER_DEG, ZERO

  ! expose integer params
  public sp,dp
  public PRINT_LEVEL, MAX_TIME_STEPS
  public NUM_SPAT_DIM, NUM_JOINTS, NUM_BODIES
  public NUM_BODY_EQN, NUM_JOINT_EQN
  public TOT_NEQN
  public NDOF_PBODY, NDOF_PJOINT     
  public TOT_NDOF_BODIES, TOT_NDOF_JOINTS, TOT_NDOF
  public ELASTIC, MBODY, CMPLX_JOINT
  public ABS_TOL, REL_TOL
  
  !=================== START OF TUNABLE PARAMS ======================!

  !-------------------------------------------------------------------!
  !  Define constants to manage precision [TUNABLE]
  !-------------------------------------------------------------------!
  integer, parameter      :: sp = kind(0.0)    ! single precision
  integer(sp), parameter  :: dp = kind(0.0d0)  ! double precision

  !-------------------------------------------------------------------!
  ! Flag used to control the printing [TUNABLE] (not functional yet)
  ! 0 = off, 1=important, 2=elaborate
  !-------------------------------------------------------------------!
  integer(sp), parameter  :: PRINT_LEVEL = 1 

  !-------------------------------------------------------------------!
  ! Specifies max time steps to steady/unsteady state [TUNABLE]
  !-------------------------------------------------------------------!
  integer(sp), parameter  :: MAX_TIME_STEPS = 100

  !-------------------------------------------------------------------!
  ! Define the system  [TUNABLE]
  !-------------------------------------------------------------------!

  ! Number of spatial dimensions (not to change yet)
  integer(sp), parameter  :: NUM_SPAT_DIM = 3

  ! Number of bodies that take part in the system (not to change yet)
  integer(sp), parameter  :: NUM_BODIES = 1

  ! Number of joints with in the system (not to change yet)
  integer(sp), parameter  :: NUM_JOINTS = 0

  !------------------------------------------------------------------!
  ! Type of simulation  [TUNABLE]
  ! (a) Rigid/ flexible 
  ! (b) Single/multibody 
  !------------------------------------------------------------------!
  ! Default: false = only rigid, only single body
  !------------------------------------------------------------------!

  ! Is it a flexible system?
  logical(sp), parameter  :: ELASTIC = .false.

  ! Is it a multi-body system?
  logical(sp), parameter  :: MBODY   = .false.

  !=================== END OF TUNABLE PARAMS =========================!

  !====================== NON-TUNABLE PARAMS =========================!
  !
  ! If these are changed the code needs to be adjusted in residual.f90,
  ! jacobian.f90, bodies.f90, joints.f90
  !
  ! The order is important too.
  !
  !-------------------------------------------------------------------!
  ! (1) NUMBER OF GOVERNING EQUATIONS
  !-------------------------------------------------------------------!

  ! 2 kinematics and 2 dynamics in vector form per body
  integer(sp), parameter  :: NUM_DYNAM_EQN = 4 

  ! Number of elastic eqns in vector form per body
  ! Note: Don't have to set this to zero for rigid case -- taken care
  ! of by the boolean ELASTIC
  integer(sp), parameter  :: NUM_ELAST_EQN = 0 

  integer(sp), parameter  :: NUM_BODY_EQN = &
       & NUM_DYNAM_EQN +  NUM_ELAST_EQN

  ! Number of joint equations in vector form per joint
  ! It this is set to zero, only single body simulation is considered
  integer(sp), parameter  :: NUM_JOINT_EQN = 0

  ! Number of governing eqns in vector form  for all bodies and joints
  integer(sp), parameter  :: TOT_NEQN = (NUM_DYNAM_EQN +&
       & NUM_ELAST_EQN)*NUM_BODIES + NUM_JOINT_EQN * NUM_JOINTS

  ! ------------------------------------------------------------------!
  !  DEGREES OF FREEDOM PER BODY AND/OR JOINT
  ! ------------------------------------------------------------------!

  ! Number of rigid DOFs per body in component form (includes dyn and 
  ! kinematic unknowns)
  integer(sp), parameter  :: NUM_RIG_DOF_PBODY = &
       & NUM_DYNAM_EQN*NUM_SPAT_DIM

  ! Number of elatic DOFs per body in component form
  integer(sp), parameter  :: NUM_ELA_DOF_PBODY = &
       & NUM_ELAST_EQN*NUM_SPAT_DIM

  ! A body has rigid and elastic degrees of freedom
  integer(sp), parameter  :: NDOF_PBODY = &
       & NUM_RIG_DOF_PBODY + NUM_ELA_DOF_PBODY

  ! Number of DOFs due to joint (usually 6)
  integer(sp), parameter  :: NDOF_PJOINT = &
       & NUM_JOINT_EQN*NUM_SPAT_DIM

  ! ------------------------------------------------------------------!
  !  DEGREES OF FREEDOM FOR 'ALL' BODIES AND/OR JOINTS
  ! ------------------------------------------------------------------!

  ! Total rigid NDOF due to all BODIES
  integer(sp), parameter  :: TOT_NUM_RIG_DOF_BODIES = &
       & NUM_RIG_DOF_PBODY*NUM_BODIES

  ! Total elastic NDOF due to all BODIES
  integer(sp), parameter  :: TOT_NUM_ELA_DOF_BODIES = &
       & NUM_ELA_DOF_PBODY*NUM_BODIES

  ! Total NDOF due to all BODIES
  integer(sp), parameter  :: TOT_NDOF_BODIES = &
       & TOT_NUM_RIG_DOF_BODIES + TOT_NUM_ELA_DOF_BODIES

  ! Total NDOF due to all JOINTS
  integer(sp), parameter  :: TOT_NDOF_JOINTS = &
       & NUM_JOINT_EQN*NUM_SPAT_DIM

  ! TOTAL number of states (DOF) from all the BODIES and JOINTS
  integer(sp), parameter  :: TOT_NDOF = &
       & TOT_NDOF_BODIES +  TOT_NDOF_JOINTS
  
  !-------------------------------------------------------------------!
  !  TYPE OF LINK (SIMPLE OR COMPLEX JOINT)
  !-------------------------------------------------------------------!
  !  (a) If joint is between two bodies, then it is a SIMPLE joint
  !  (b) If joint is between >two bodies, then complex joints exist
  !
  ! Note: Why is this important?
  !
  ! This will affect the sparsity of the linear system and we would
  ! like to maintain sparsity by appropriate ordering of the entities
  !-------------------------------------------------------------------!
 
  logical, parameter :: CMPLX_JOINT = .false.   
  
  !================= END OF NON-TUNABLE PARAMS =======================!

  !-------------------------------------------------------------------!
  ! Define double constants [TUNABLE]
  !-------------------------------------------------------------------!
  real(dp), parameter     :: PI = 22.0_dp/7.0_dp
  real(dp), parameter     :: DEG_PER_RAD = 180.0_dp/PI
  real(dp), parameter     :: RAD_PER_DEG = PI/180.0_dp
  real(dp), parameter     :: ZERO = 0.0_dp
  
  ! tolerances used in time-integration
  real(dp), parameter     :: REL_TOL = 1.e-6_dp
  real(dp), parameter     :: ABS_TOL = 1.e-6_dp

end module global_constants
