!=====================================================================!
! Rigid body implementation of the abstract body class
!---------------------------------------------------------------------!
! Has functions to:
! (a) create rigid body based on supplied state, mass, J, c etc
! (b) update the rigid body state variables alone during time-marching
! (d) compute the rotation and angular rate matrices of the body
! (c) 'toString' like function to print the body props (state+attrs)
!=====================================================================!

module rigid_body_class

  ! module references
  use body_class
  use global_constants
  use global_variables
  use types, only: matrix, vector
  use utils
  use dispmodule, only : disp

  ! module options
  implicit none
  private
  public :: rigid_body, print_rigid_body
  public :: get_rotation
  public :: get_angrate, get_angrate_dot, get_angrate_inv

  !*******************************************************************!
  ! RIGID_BODY datatype can be used to fully characterize the STATE and 
  ! ATTRIBUTES of a rigid body . A rigid_body object contains virtually
  ! everything about the body 
  !*******************************************************************!
  
  type, extends(body) :: rigid_body

     !----------------------------------------------------------------!
     ! rigid body state variables
     !----------------------------------------------------------------!

     ! origin of the body frame
     type(vector) :: r              

     ! orientation of the body frame with inertial (euler angles)
     type(vector) :: theta          

     ! velocity of the origin with respect to inertial
     type(vector) :: v              

     ! angular velocity of the body frame with respect to inertial
     type(vector) :: omega          

     !----------------------------------------------------------------!
     ! time derivative of states
     !----------------------------------------------------------------!
     type(vector) :: r_dot
     type(vector) :: theta_dot
     type(vector) :: v_dot
     type(vector) :: omega_dot
     !----------------------------------------------------------------!
     ! Body Attributes
     !----------------------------------------------------------------!

     real(dp)     :: mass    ! mass (denoted m in paper)   

     !     The format for c is: (in body frame)
     !     c = [ cx,  cy,  cz ]

     type(vector) :: c              ! first moment of inertia

     !  The format for J is: (in body frame)
     !  J = [ Jxx,  Jxy,  Jxz ] = [ J[0],  J[1],  J[2] ]
     !  . = [    ,  Jyy,  Jyz ] = [     ,  J[3],  J[4] ]
     !  . = [    ,     ,  Jzz ] = [     ,      ,  J[5] ]

     type(matrix) :: J              ! second moment of inertia

     type(matrix) :: C_mat          ! rotation matrix
     type(matrix) :: S
     type(matrix) :: S_dot          ! transformation matrix

     type(vector) :: fr             ! reaction force
     type(vector) :: gr             ! reaction torque

     type(vector) :: rj             ! position of joint
     type(vector) :: g              ! gravity vector in local frame

     real(dp)     :: KE             ! kinetic energy of the body
     real(dp)     :: PE             ! potential energy of the body

  end type rigid_body

!*******************************************************************!
! A common interface for creating bodies and updating existing ones
!-------------------------------------------------------------------!
! The inputs are described in the methods under the interface. One 
! can use the same interface for updating existing bodies.
! Shortly: use this to create and update bodies
!*******************************************************************!

interface rigid_body

  procedure constructor

end interface rigid_body

!*******************************************************************!
! A common interface for different ways of getting rotation matrix
!-------------------------------------------------------------------!
! (a) get_rotation_from_angles_vec   -- > input theta as VECTOR
! (b) get_rotation_from_angles_array -- > input theta(3) as array
! (c) get_rotation_from_cosines -- > input dir cosines and sines
!*******************************************************************!
interface get_rotation
  module procedure get_rotation_from_angles_array, &
       &get_rotation_from_angles_vec, get_rotation_from_cosines
end interface get_rotation

!*******************************************************************!
! A common interface for different ways of getting ang rate  matrix
!-------------------------------------------------------------------!
! (a) get_angrate_from_angles_vec    -- > input VECTOR theta
! (b) get_angrate_from_angles_array  -- > input theta(3)
! (c) get_angrate_from_cosines       -- > input dir cosines and sines
!*******************************************************************!
interface get_angrate
  module procedure get_angrate_from_angles_vec, &
       &get_angrate_from_angles_array, get_angrate_from_cosines
end interface get_angrate

!*******************************************************************!
! A common interface for different ways of getting the time 
! derivative of the ang rate matrix
!-------------------------------------------------------------------!
! (a) get_angrate_dot_vec -- > input VECTOR theta, theta_dot
! (b) get_angrate_dot_array -- > input theta(3), theta_dot(3)
! (c) get_angrate_dot_cosines  -- > input dir cosines and sines
!*******************************************************************!
interface get_angrate_dot
  module procedure  get_angrate_dot_array, &
       &get_angrate_dot_vec, get_angrate_dot_cosines
end interface get_angrate_dot

!*******************************************************************!
! A common interface for different ways of getting the inverse
! of the ang rate matrix
!-------------------------------------------------------------------!
! (a) get_angrate_inv_vec -- > input VECTOR theta, theta_dot
! (b) get_angrate_inv_array -- > input theta(3), theta_dot(3)
! (c) get_angrate_inv_cosines  -- > input dir cosines and sines
!*******************************************************************!
interface get_angrate_inv
  module procedure get_angrate_inv_vec, &
       &get_angrate_inv_array, &
       &get_angrate_inv_cosines
end interface get_angrate_inv

contains

  !*******************************************************************!
  ! create a body using the supplied parameters
  !-------------------------------------------------------------------!
  ! mass   : mass of the body
  ! c      : first moment of mass 
  ! J      : second moment of mass
  ! fr     : reaction force
  ! gr     : reaction torque
  ! q, qdot: state vector and time derivatives
  ! qddot  : second time derivative of the state (used in elastic only)
  !*******************************************************************!
  
  function constructor(mass, c, J, fr, gr, q, qdot) result(this)

    ! inputs
    real(dp), intent(in) :: mass 
    real(dp), intent(in) :: c(NUM_SPAT_DIM)
    real(dp), intent(in) :: J(NUM_SPAT_DIM, NUM_SPAT_DIM)
    real(dp), intent(in) :: fr(NUM_SPAT_DIM)
    real(dp), intent(in) :: gr(NUM_SPAT_DIM)

    real(dp), intent(in) :: q(NDOF_PBODY)
    real(dp), intent(in) :: qdot(NDOF_PBODY)

    ! input/output
    type(rigid_body) :: this

    !-----------------------------------------------------------------!
    ! Inertial properties of the body
    !-----------------------------------------------------------------!

    ! mass of the body
    this%mass = mass

    ! moment of inertia in body-fixed frame
    this%J = matrix(J)       !-mass*skew(re)*skew(re)

    !first moment of inertia in body-fixed frame: mass*(cg location)
    this%c = vector(c)       ! mass*re

    !-----------------------------------------------------------------!
    ! set the state into the body
    !-----------------------------------------------------------------!

    this%r         = vector(q(1:3))
    this%theta     = vector(q(4:6))
    this%v         = vector(q(7:9))
    this%omega     = vector(q(10:12))

    !-----------------------------------------------------------------!
    ! set the time derivatives of state into the body
    !-----------------------------------------------------------------!


    this%r_dot     = vector(qdot(1:3))
    this%theta_dot = vector(qdot(4:6))
    this%v_dot     = vector(qdot(7:9))
    this%omega_dot = vector(qdot(10:12))


    !-----------------------------------------------------------------!
    ! update the rotation and angular rate matrices
    !-----------------------------------------------------------------!

    this%C_mat     = get_rotation(this%theta)
    this%S         = get_angrate(this%theta)
    this%S_dot     = get_angrate_dot(this%theta, this%theta_dot)

    !-----------------------------------------------------------------!
    ! update the new direction of the gravity vector in body frame
    !-----------------------------------------------------------------!

    this%g = this%C_mat*GRAV

    !-----------------------------------------------------------------!
    ! Mechanical Energy 
    !-----------------------------------------------------------------!

    ! update the kinetic energy
    this%KE = 0.5_dp*(this%mass * this%v *  this%v &
         & + this%omega*this%J* this%omega) ! + coupling term

    ! update potential energy
    this%PE = this%mass*this%g*this%r

    !-----------------------------------------------------------------!
    ! Joint reactions
    !-----------------------------------------------------------------!

    ! reaction force 
    this%fr    = vector(fr)   

    ! reaction torque
    this%gr    = vector(gr)

    if (this%mass .eq. ZERO) stop "ERROR: Body with zero mass??!?"

    !call print_body(this)

  end function constructor
  
  !*******************************************************************!
  ! Returns the rotation matrix based on the euler angles
  ! Compute the 3-2-1 Euler angle rotation
  
  ! C = C1(theta_1)*C2(theta_2)*C3(theta_3)
  !
  ! Input: theta as an array
  ! Output: CMAT of type MATRIX
  ! 
  ! Ref: Section 2.2 Eq. 18/19 Hughes
  !*******************************************************************!
  function get_rotation_from_angles_array(theta) result(CMAT)

    real(dp), intent(in)       :: theta(NUM_SPAT_DIM)    
    type(matrix)               :: CMAT
    real(dp)                   :: c1, c2, c3, s1, s2, s3

    ! Compute the sin/cos of the Euler angles
    c1 = cos(theta(1))
    s1 = sin(theta(1))

    c2 = cos(theta(2))
    s2 = sin(theta(2))

    c3 = cos(theta(3))
    s3 = sin(theta(3))

    CMAT = get_rotation_from_cosines(c1, c2, c3, s1, s2, s3)

  end function get_rotation_from_angles_array

  !*******************************************************************!
  ! Returns the rotation matrix based on the euler angles
  !
  ! Compute the 3-2-1 Euler angle rotation
  !
  ! CMAT = C1(theta_1)*C2(theta_2)*C3(theta_3)
  !
  ! Input: theta of type VECTOR
  ! Output: CMAT of type MATRIX
  ! 
  ! Ref: Section 2.2 Eq. 18/19 Hughes
  !*******************************************************************!
  function get_rotation_from_angles_vec(thetain) result(CMAT)

    type(vector), intent(in)   :: thetain
    real(dp)                   :: theta(NUM_SPAT_DIM)
    type(matrix)               :: CMAT

    ! covert to array form
    theta = array(thetain)

    ! call the method that takes angles array
    CMAT  =  get_rotation_from_angles_array(theta)

  end function get_rotation_from_angles_vec

  !*******************************************************************!
  ! Returns the rotation matrix (euler angles) based on the sines and 
  ! cosines of the euler angles
  !
  ! Compute the 3-2-1 Euler angle rotation

  ! Input: sines and cosines of the euler angles
  ! Output: CMAT of type MATRIX
  !
  ! Ref: Section 2.2 Eq. 18/19 Hughes
  !*******************************************************************!
  function get_rotation_from_cosines(c1,c2,c3,s1,s2,s3) result(CMAT)

    real(dp), intent(in)       :: c1, c2, c3, s1, s2, s3
    type(matrix)               :: CMAT

    CMAT = matrix((/ c2*c3, c2*s3, -s2,&
         & s1*s2*c3 - c1*s3, s1*s2*s3 + c1*c3, s1*c2,&
         & c1*s2*c3 + s1*s3, c1*s2*s3 - s1*c3, c1*c2 /))

  end function get_rotation_from_cosines


  !*******************************************************************!
  ! Returns the ang rate matrix from the supplied euler angles
  !
  ! Compute the 3-2-1 Euler angle rotation. The S matrix does not 
  ! depend on c3/s3, however we keep these as inputs in case we ever
  ! want to change the Euler sequence.
  !
  ! Input : theta of type VECTOR
  ! Output: SMAT  of type MATRIX
  ! 
  ! Ref: Section 2.3 Eq. 24/25 Hughes
  ! ******************************************************************!
  function get_angrate_from_angles_vec(thetain) result(SMAT)

    type(vector) :: thetain
    real(dp)     :: theta(NUM_SPAT_DIM)    
    type(matrix) :: SMAT

    ! convert to array
    theta = array(thetain)

    ! call the function with array signature
    SMAT = get_angrate_from_angles_array(theta)

  end function get_angrate_from_angles_vec

  !*******************************************************************!
  ! Returns the ang rate matrix from the supplied euler angles
  !
  ! Compute the 3-2-1 Euler angle rotation. The S matrix does not 
  ! depend on c3/s3, however we keep these as inputs in case we ever
  ! want to change the Euler sequence.
  !
  ! Input : theta as an array
  ! Output: SMAT  of type MATRIX
  ! 
  ! Ref: Section 2.3 Eq. 24/25 Hughes
  ! ******************************************************************!
  function get_angrate_from_angles_array(theta) result(SMAT)

    real(dp)     :: theta(NUM_SPAT_DIM)    
    type(matrix) :: SMAT
    real(dp)     :: c1, c2, c3, s1, s2, s3

    ! Compute the sin/cos of the Euler angles
    c1 = cos(theta(1))
    s1 = sin(theta(1))

    c2 = cos(theta(2))
    s2 = sin(theta(2))

    c3 = cos(theta(3))
    s3 = sin(theta(3))

    SMAT = get_angrate_from_cosines(c1, c2, c3, s1, s2, s3)

  end function get_angrate_from_angles_array

  !*******************************************************************!
  ! Returns the rotation matrix (euler angles) based on the sines and 
  ! cosines of the euler angles
  !
  ! Compute the 3-2-1 Euler angle rotation. The S matrix does not 
  ! depend on c3/s3, however we keep these as inputs in case we ever
  ! want to change the Euler sequence.
  !
  ! Input : sines and cosines of euler angles
  ! Output: SMAT of type MATRIX
  ! 
  ! Ref: Section 2.3 Eq. 24/25 Hughes
  !*******************************************************************!
  function get_angrate_from_cosines(c1,c2,c3,s1,s2,s3) result(SMAT)

    real(dp), intent(in)       :: c1, c2, c3, s1, s2, s3
    type(matrix)               :: SMAT

    SMAT = matrix((/ 1.0_dp, 0.0_dp, -s2, &
         & 0.0_dp,  c1,  s1*c2, &
         & 0.0_dp,  -s1,  c1*c2 /))

  end function get_angrate_from_cosines

  !-----------------------------------------------------------
  ! Returns the time derivative of angular rate matrix
  !
  ! we use the 3-2-1 Euler angles, the S matrix does not depend
  ! on c3/s3, however we keep these as inputs in case we ever  
  ! want to change the Euler sequence.
  !
  ! Input : the euler angles and time derivatives as VECTOR
  ! Output: SMAT_DOT
  !
  !-----------------------------------------------------------
  function get_angrate_dot_vec(thetain, dthetain) result(SMAT_DOT)

    type(vector) :: thetain, dthetain
    type(matrix) :: SMAT_DOT

    real(dp)     :: theta(NUM_SPAT_DIM), dtheta(NUM_SPAT_DIM)

    ! convert vec to array
    theta = array(thetain); dtheta=array(dthetain);   

    ! call the function matching array signature
    SMAT_DOT = get_angrate_dot_array(theta,dtheta)

  end function get_angrate_dot_vec


  !-----------------------------------------------------------
  ! Returns the time derivative of angular rate matrix
  !
  ! we use the 3-2-1 Euler angles, the S matrix does not depend
  ! on c3/s3, however we keep these as inputs in case we ever  
  ! want to change the Euler sequence.
  !
  ! Input : the euler angles and time derivatives as arrays
  ! Output: SMAT_DOT
  !
  !-----------------------------------------------------------
  function get_angrate_dot_array(theta, dtheta) result(SMAT_DOT)

    real(dp)     :: theta(NUM_SPAT_DIM), dtheta(NUM_SPAT_DIM)
    type(matrix) :: SMAT_DOT
    real(dp)     :: c1, c2, c3, s1, s2, s3

    ! Compute the sin/cos of the Euler angles

    c1 = cos(theta(1))
    s1 = sin(theta(1))

    c2 = cos(theta(2))
    s2 = sin(theta(2))

    c3 = cos(theta(3))
    s3 = sin(theta(3))

    SMAT_DOT = get_angrate_dot_cosines( c1, c2, c3, s1, s2, s3, dtheta)

  end function get_angrate_dot_array


  !-----------------------------------------------------------
  ! Returns the time derivative of angular rate matrix
  !
  ! we use the 3-2-1 Euler angles, the S matrix does not depend
  ! on c3/s3, however we keep these as inputs in case we ever  
  ! want to change the Euler sequence.
  !
  ! Input : The sines and cosines of euler angles
  ! Output: SMAT_DOT
  !
  !-----------------------------------------------------------
  function get_angrate_dot_cosines( c1, c2, c3, s1, s2, s3, dtheta) &
       &result(SMAT_DOT)

    real(dp)     :: dtheta(NUM_SPAT_DIM) 
    real(dp)     :: c1, c2, c3, s1, s2, s3
    type(matrix) :: SMAT_DOT

    SMAT_DOT= matrix( (/ 0.0_dp,  0.0_dp,  -c2*dtheta(2), &
         & 0.0_dp, -s1*dtheta(1), c1*c2*dtheta(1)-s1*s2*dtheta(2),&
         & 0.0_dp, -c1*dtheta(1), -s1*c2*dtheta(1)-c1*s2*dtheta(2)/))

  end function get_angrate_dot_cosines


  ! ******************************************************************!
  ! Returns the inverse of the angular rate matrix for the supplied
  ! theta vector
  ! ******************************************************************!
  function get_angrate_inv_vec(thetain) &
       &result(SMAT_INV)

    type(vector) :: thetain
    type(matrix) :: SMAT_INV
    real(dp)     :: theta(NUM_SPAT_DIM)

    ! decompose the vector into array
    theta = array(thetain)

    ! call the method that takes array as input
    SMAT_INV =  get_angrate_inv_array(theta)

  end function get_angrate_inv_vec

  ! ******************************************************************!
  ! Returns the inverse of the angular rate matrix for the supplied
  ! theta array
  ! ******************************************************************!
  function get_angrate_inv_array(theta) &
       &result(SMAT_INV)

    real(dp)     :: theta(NUM_SPAT_DIM)
    real(dp)     :: c1, c2, c3, s1, s2, s3
    type(matrix) :: SMAT_INV

    c1 = cos(theta(1))
    s1 = sin(theta(1))

    c2 = cos(theta(2))
    s2 = sin(theta(2))

    c3 = cos(theta(3))
    s3 = sin(theta(3))

    SMAT_INV =  get_angrate_inv_cosines(  c1, c2, c3, s1, s2, s3)

  end function get_angrate_inv_array

  ! ******************************************************************!
  ! Returns the inverse of the angular rate matrix for the supplied
  ! direction cosines
  ! ******************************************************************!
  function get_angrate_inv_cosines( c1, c2, c3, s1, s2, s3) &
       &result(SMAT_INV)

    real(dp)     :: c1, c2, c3, s1, s2, s3
    type(matrix) :: SMAT_INV

    SMAT_INV= matrix( (/  1.0_dp, s1*s2/c2, c1*s2/c2, &
         &0.0_dp, c1, -s1, 0.0_dp,&
         & s1/c2, c1/c2 /))

  end function get_angrate_inv_cosines

  !*******************************************************************!
  ! routine that prints the state and properties of the body
  !******************************************************************!
  subroutine print_rigid_body(this)

    class(rigid_body):: this

    call disp('======================================================')
    call disp('---------------------BODY----------------------------')
    call disp('======================================================')
    call disp('')
    call disp('   > r          =   ', array(this%r), SEP=', ', &
         &ORIENT = 'ROW')
    call disp('   > theta      =   ', array(this%theta), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('   > v          =   ', array(this%v), SEP=', ', &
         &ORIENT = 'ROW')
    call disp('   > omega      =   ', array(this%omega), SEP=', ',&
         & ORIENT = 'ROW')
    call DISP('')
    call disp('   > r_dot      =   ', array(this%r_dot), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('   > theta_dot  =   ', array(this%theta_dot), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('   > v_dot      =   ', array(this%v_dot), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('   > omega_dot  =   ', array(this%omega_dot), SEP=', ',&
         &ORIENT = 'ROW')
    call DISP('')
    call disp('   > c          =   ', array(this%c), SEP=', ',&
         & ORIENT = 'ROW')
    call DISP('')
    call disp('   > J          =   ', matrix(this%J))
    call DISP('')
    call disp('   > Rot Mat    =   ', matrix(this%C_mat))
    call DISP('')
    call disp('   > Angrt Mat  =   ', matrix(this%S))
    call DISP('')
    call disp('   > S_dot      =   ', matrix(this%S_dot))
    call DISP('')
    call disp('   > fr         =   ', array(this%fr), SEP=', ', &
         &ORIENT = 'ROW')

    call disp('   > gr         =   ', array(this%gr), SEP=', ', &
         &ORIENT = 'ROW')
    call disp('')
    call disp('   > gravity    =   ', array(this%g), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('')
    call disp('   > Pot Energy =   ', this%PE)
    call disp('   > Kin Energy =   ', this%KE)
    call disp('   > Tot Energy =   ', this%KE + this%PE)
    call disp('======================================================')
  end subroutine print_rigid_body

end module rigid_body_class
