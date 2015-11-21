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

  !  private

  !  public :: rigid_body , print_rigid_body

  !  public :: find_rotation
  !  public :: find_angrate, find_angrate_dot, find_angrate_inv

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
     !  J = [ Jxx,  Jxy,  Jxz ] = [ J[1],  J[2],  J[3] ]
     !  . = [    ,  Jyy,  Jyz ] = [     ,  J[4],  J[5] ]
     !  . = [    ,     ,  Jzz ] = [     ,      ,  J[6] ]

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

   contains

     ! type-bound getters
     procedure:: get_r, get_theta, get_v, get_omega
     procedure:: get_r_dot, get_theta_dot, get_v_dot, get_omega_dot
     procedure:: get_mass, get_first_moment, get_second_moment
     procedure:: get_reaction_force, get_reaction_torque
     procedure:: get_potential_energy, get_kinetic_energy
     procedure:: get_joint_location, get_gravity
     procedure:: get_rotation, get_angrate, get_angrate_dot

     ! type-bound setters
     procedure:: set_r, set_theta, set_v, set_omega
     procedure:: set_r_dot, set_theta_dot, set_v_dot, set_omega_dot
     procedure:: set_mass, set_first_moment, set_second_moment
     procedure:: set_reaction_force, set_reaction_torque
     procedure:: set_potential_energy, set_kinetic_energy
     procedure:: set_joint_location, set_gravity
     procedure:: set_rotation, set_angrate, set_angrate_dot
     


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
  ! (a) find_rotation_from_angles_vec   -- > input theta as VECTOR
  ! (b) find_rotation_from_angles_array -- > input theta(3) as array
  ! (c) find_rotation_from_cosines -- > input dir cosines and sines
  !*******************************************************************!

  interface find_rotation
     module procedure &
          & find_rotation_from_angles_array, &
          & find_rotation_from_angles_vec,&
          & find_rotation_from_cosines
  end interface find_rotation

  !*******************************************************************!
  ! A common interface for different ways of getting ang rate  matrix
  !-------------------------------------------------------------------!
  ! (a) find_angrate_from_angles_vec    -- > input VECTOR theta
  ! (b) find_angrate_from_angles_array  -- > input theta(3)
  ! (c) find_angrate_from_cosines       -- > input dir cosines and sines
  !*******************************************************************!

  interface find_angrate
     module procedure &
          & find_angrate_from_angles_vec, &
          & find_angrate_from_angles_array,&
          & find_angrate_from_cosines
  end interface find_angrate

  !*******************************************************************!
  ! A common interface for different ways of getting the time 
  ! derivative of the ang rate matrix
  !-------------------------------------------------------------------!
  ! (a) find_angrate_dot_vec -- > input VECTOR theta, theta_dot
  ! (b) find_angrate_dot_array -- > input theta(3), theta_dot(3)
  ! (c) find_angrate_dot_cosines  -- > input dir cosines and sines
  !*******************************************************************!

  interface find_angrate_dot
     module procedure &
          & find_angrate_dot_array, &
          & find_angrate_dot_vec, &
          & find_angrate_dot_cosines
  end interface find_angrate_dot

  !*******************************************************************!
  ! A common interface for different ways of getting the inverse
  ! of the ang rate matrix
  !-------------------------------------------------------------------!
  ! (a) find_angrate_inv_vec -- > input VECTOR theta, theta_dot
  ! (b) find_angrate_inv_array -- > input theta(3), theta_dot(3)
  ! (c) find_angrate_inv_cosines  -- > input dir cosines and sines
  !*******************************************************************!

  interface find_angrate_inv
     module procedure &
          & find_angrate_inv_vec, &
          & find_angrate_inv_array, &
          & find_angrate_inv_cosines
  end interface find_angrate_inv

contains

  !*******************************************************************!
  ! Getter for the position of the body
  !*******************************************************************!

  function get_r(this)

    class(rigid_body) :: this
    type(vector)     :: get_r

    get_r = this % r

  end function get_r

  !*******************************************************************!
  ! Getter for the orientation of the body
  !*******************************************************************!

  function get_theta(this)

    class(rigid_body) :: this
    type(vector)     :: get_theta

    get_theta = this % theta

  end function get_theta

  !*******************************************************************!
  ! Getter for the velocity of the body
  !*******************************************************************!

  function get_v(this)

    class(rigid_body) :: this
    type(vector)     :: get_v

    get_v = this % v

  end function get_v

  !*******************************************************************!
  ! Getter for the angular velocity of the body
  !*******************************************************************!

  function get_omega(this)

    class(rigid_body) :: this
    type(vector)     :: get_omega

    get_omega = this % omega

  end function get_omega

  !*******************************************************************!
  ! Getter for the kinetic energy of the body
  !*******************************************************************!

  function get_kinetic_energy(this)

    class(rigid_body) :: this
    real(dp) :: get_kinetic_energy

    get_kinetic_energy = this % KE

  end function get_kinetic_energy

  !*******************************************************************!
  ! Getter for the potential energy of the body
  !*******************************************************************!

  function get_potential_energy(this)

    class(rigid_body) :: this
    real(dp)     :: get_potential_energy

    get_potential_energy = this % PE

  end function get_potential_energy

  !*******************************************************************!
  ! Getter for the rdot of the body
  !*******************************************************************!

  function get_r_dot(this)

    class(rigid_body) :: this
    type(vector)     :: get_r_dot

    get_r_dot = this % r_dot

  end function get_r_dot

  !*******************************************************************!
  ! Getter for theta_dot of the body
  !*******************************************************************!

  function get_theta_dot(this)

    class(rigid_body) :: this
    type(vector)     :: get_theta_dot

    get_theta_dot = this % theta_dot

  end function get_theta_dot

  !*******************************************************************!
  ! Getter for vdot of the body
  !*******************************************************************!

  function get_v_dot(this)

    class(rigid_body) :: this
    type(vector)     :: get_v_dot

    get_v_dot = this % v_dot

  end function get_v_dot

  !*******************************************************************!
  ! Getter for the angular velocity of the body
  !*******************************************************************!

  function get_omega_dot(this)

    class(rigid_body) :: this
    type(vector)     :: get_omega_dot

    get_omega_dot = this % omega_dot

  end function get_omega_dot


  !*******************************************************************!
  ! Getter for the mass of the body
  !*******************************************************************!

  function get_mass(this)

    class(rigid_body) :: this
    real(dp)     :: get_mass

    get_mass = this % mass

  end function get_mass

  !*******************************************************************!
  ! Getter for the first_moment of the body
  !*******************************************************************!

  function get_first_moment (this)

    class(rigid_body) :: this
    type(vector)      :: get_first_moment

    get_first_moment = this % c

  end function get_first_moment

  !*******************************************************************!
  ! Getter for the second_moment of the body
  !*******************************************************************!

  function get_second_moment (this)

    class(rigid_body) :: this
    type(matrix)      :: get_second_moment

    get_second_moment = this % J

  end function get_second_moment

  !*******************************************************************!
  ! Getter for the reaction_force of the body
  !*******************************************************************!

  function get_reaction_force(this)

    class(rigid_body) :: this
    type(vector) :: get_reaction_force

    get_reaction_force = this % fr

  end function get_reaction_force
  

  !*******************************************************************!
  ! Getter for the reaction_torque of the body
  !*******************************************************************!

  function get_reaction_torque(this)

    class(rigid_body) :: this
    type(vector) :: get_reaction_torque

    get_reaction_torque = this % gr

  end function get_reaction_torque


  !*******************************************************************!
  ! Getter for the joint_location on the body
  !*******************************************************************!
  
  function get_joint_location(this)

    class(rigid_body) :: this
    type(vector) :: get_joint_location

    get_joint_location = this % rj

  end function get_joint_location

  !*******************************************************************!
  ! Getter for the gravity vector in the body
  !*******************************************************************!
  
  function get_gravity(this)

    class(rigid_body) :: this
    type(vector) :: get_gravity

    get_gravity = this % g

  end function get_gravity

  !*******************************************************************!
  ! Getter for body to inertial rotation matrix
  !*******************************************************************!

  function get_rotation(this)

    class(rigid_body) :: this
    type(matrix) :: get_rotation

    get_rotation = this % C_mat
    
  end function get_rotation


  !*******************************************************************!
  ! Getter for the angular rate matrix
  !*******************************************************************!

  function get_angrate(this)

    class(rigid_body) :: this
    type(matrix) :: get_angrate

    get_angrate = this % S

  end function get_angrate


  !*******************************************************************!
  ! Getter for the angrate_dot matrix
  !*******************************************************************!

  function get_angrate_dot(this)

    class(rigid_body) :: this
    type(matrix) :: get_angrate_dot

    get_angrate_dot = this % S_dot

  end function get_angrate_dot
  
  !===================== SETTERS TO TYPE VARIABLES ===================!
 
  !*******************************************************************!
  ! Setter for the position of the body
  !*******************************************************************!

  subroutine set_r(this, r)

    class(rigid_body) :: this
    type(vector)     :: r

    this % r = r

  end subroutine set_r

  !*******************************************************************!
  ! Setter for the orientation of the body
  !*******************************************************************!

  subroutine set_theta(this, theta)

    class(rigid_body) :: this
    type(vector)     :: theta

    this % theta = theta

  end subroutine set_theta

  !*******************************************************************!
  ! Setter for the velocity of the body
  !*******************************************************************!

  subroutine set_v(this, v)

    class(rigid_body) :: this
    type(vector)     :: v

    this % v = v

  end subroutine set_v

  !*******************************************************************!
  ! Setter for the angular velocity of the body
  !*******************************************************************!

  subroutine set_omega(this, omega)

    class(rigid_body) :: this
    type(vector)      :: omega

    this % omega = omega

  end subroutine set_omega

  !*******************************************************************!
  ! Setter for the kinetic energy of the body
  !*******************************************************************!

  subroutine set_kinetic_energy(this, KE)

    class(rigid_body) :: this
    real(dp) :: KE

    this % KE = KE

  end subroutine set_kinetic_energy

  !*******************************************************************!
  ! Setter for the potential energy of the body
  !*******************************************************************!

  subroutine set_potential_energy(this, PE)

    class(rigid_body) :: this
    real(dp)     :: PE

    this % PE = PE

  end subroutine set_potential_energy

  !*******************************************************************!
  ! Setter for the rdot of the body
  !*******************************************************************!

  subroutine set_r_dot(this, r_dot)

    class(rigid_body) :: this
    type(vector)      :: r_dot

    this % r_dot = r_dot

  end subroutine set_r_dot

  !*******************************************************************!
  ! Setter for theta_dot of the body
  !*******************************************************************!

  subroutine set_theta_dot(this, theta_dot)

    class(rigid_body) :: this
    type(vector)      :: theta_dot

    this % theta_dot = theta_dot

  end subroutine set_theta_dot

  !*******************************************************************!
  ! Setter for vdot of the body
  !*******************************************************************!

  subroutine set_v_dot(this, v_dot)

    class(rigid_body) :: this
    type(vector)      :: v_dot

    this % v_dot = v_dot

  end subroutine set_v_dot

  !*******************************************************************!
  ! Setter for the angular velocity of the body
  !*******************************************************************!

  subroutine set_omega_dot(this, omega_dot)

    class(rigid_body) :: this
    type(vector)      :: omega_dot

    this % omega_dot = omega_dot

  end subroutine set_omega_dot


  !*******************************************************************!
  ! Setter for the mass of the body
  !*******************************************************************!

  subroutine set_mass(this, mass)

    class(rigid_body) :: this
    real(dp)      :: mass

    this % mass = mass

  end subroutine set_mass

  !*******************************************************************!
  ! Setter for the first_moment of the body
  !*******************************************************************!

  subroutine set_first_moment (this, c)

    class(rigid_body) :: this
    type(vector)      :: c

    this % c = c

  end subroutine set_first_moment

  !*******************************************************************!
  ! Setter for the second_moment of the body
  !*******************************************************************!

  subroutine set_second_moment (this, J)

    class(rigid_body) :: this
    type(matrix)      :: J

    this % J = J

  end subroutine set_second_moment
  
  !*******************************************************************!
  ! Setter for the reaction_force of the body
  !*******************************************************************!

  subroutine set_reaction_force(this, fr)

    class(rigid_body) :: this
    type(vector) :: fr

    this % fr = fr

  end subroutine set_reaction_force
  

  !*******************************************************************!
  ! Setter for the reaction_torque of the body
  !*******************************************************************!

  subroutine set_reaction_torque(this, gr)

    class(rigid_body) :: this
    type(vector) :: gr

    this % gr = gr

  end subroutine set_reaction_torque


  !*******************************************************************!
  ! Setter for the joint_location on the body
  !*******************************************************************!
  
  subroutine set_joint_location(this, rj)

    class(rigid_body) :: this
    type(vector) :: rj

    this % rj = rj

  end subroutine set_joint_location

  !*******************************************************************!
  ! Setter for the gravity vector in the body
  !*******************************************************************!
  
  subroutine set_gravity(this, g)

    class(rigid_body) :: this
    type(vector) :: g

    this % g = g

  end subroutine set_gravity

  !*******************************************************************!
  ! Setter for body to inertial rotation matrix
  !*******************************************************************!

  subroutine set_rotation(this, C_mat)

    class(rigid_body) :: this
    type(matrix) :: C_mat

    this % C_mat = C_mat
    
  end subroutine set_rotation

  !*******************************************************************!
  ! Setter for the angular rate matrix
  !*******************************************************************!
  
  subroutine set_angrate(this, S)

    class(rigid_body) :: this
    type(matrix) :: S

    this % S = S

  end subroutine set_angrate

  !*******************************************************************!
  ! Setter for the angrate_dot matrix
  !*******************************************************************!

  subroutine set_angrate_dot(this, S_dot)

    class(rigid_body) :: this
    type(matrix) :: S_dot

    this % S_dot = S_dot

  end subroutine set_angrate_dot

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

    this%C_mat     = find_rotation(this%theta)
    this%S         = find_angrate(this%theta)
    this%S_dot     = find_angrate_dot(this%theta, this%theta_dot)

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
  function find_rotation_from_angles_array(theta) result(CMAT)

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

    CMAT = find_rotation_from_cosines(c1, c2, c3, s1, s2, s3)

  end function find_rotation_from_angles_array

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
  function find_rotation_from_angles_vec(thetain) result(CMAT)

    type(vector), intent(in)   :: thetain
    real(dp)                   :: theta(NUM_SPAT_DIM)
    type(matrix)               :: CMAT

    ! covert to array form
    theta = array(thetain)

    ! call the method that takes angles array
    CMAT  =  find_rotation_from_angles_array(theta)

  end function find_rotation_from_angles_vec

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
  function find_rotation_from_cosines(c1,c2,c3,s1,s2,s3) result(CMAT)

    real(dp), intent(in)       :: c1, c2, c3, s1, s2, s3
    type(matrix)               :: CMAT

    CMAT = matrix((/ c2*c3, c2*s3, -s2,&
         & s1*s2*c3 - c1*s3, s1*s2*s3 + c1*c3, s1*c2,&
         & c1*s2*c3 + s1*s3, c1*s2*s3 - s1*c3, c1*c2 /))

  end function find_rotation_from_cosines


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
  function find_angrate_from_angles_vec(thetain) result(SMAT)

    type(vector) :: thetain
    real(dp)     :: theta(NUM_SPAT_DIM)    
    type(matrix) :: SMAT

    ! convert to array
    theta = array(thetain)

    ! call the function with array signature
    SMAT = find_angrate_from_angles_array(theta)

  end function find_angrate_from_angles_vec

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
  function find_angrate_from_angles_array(theta) result(SMAT)

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

    SMAT = find_angrate_from_cosines(c1, c2, c3, s1, s2, s3)

  end function find_angrate_from_angles_array

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
  function find_angrate_from_cosines(c1,c2,c3,s1,s2,s3) result(SMAT)

    real(dp), intent(in)       :: c1, c2, c3, s1, s2, s3
    type(matrix)               :: SMAT

    SMAT = matrix((/ 1.0_dp, 0.0_dp, -s2, &
         & 0.0_dp,  c1,  s1*c2, &
         & 0.0_dp,  -s1,  c1*c2 /))

  end function find_angrate_from_cosines

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
  function find_angrate_dot_vec(thetain, dthetain) result(SMAT_DOT)

    type(vector) :: thetain, dthetain
    type(matrix) :: SMAT_DOT

    real(dp)     :: theta(NUM_SPAT_DIM), dtheta(NUM_SPAT_DIM)

    ! convert vec to array
    theta = array(thetain); dtheta=array(dthetain);   

    ! call the function matching array signature
    SMAT_DOT = find_angrate_dot_array(theta,dtheta)

  end function find_angrate_dot_vec


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
  function find_angrate_dot_array(theta, dtheta) result(SMAT_DOT)

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

    SMAT_DOT = find_angrate_dot_cosines( c1, c2, c3, s1, s2, s3, dtheta)

  end function find_angrate_dot_array


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
  function find_angrate_dot_cosines( c1, c2, c3, s1, s2, s3, dtheta) &
       &result(SMAT_DOT)

    real(dp)     :: dtheta(NUM_SPAT_DIM) 
    real(dp)     :: c1, c2, c3, s1, s2, s3
    type(matrix) :: SMAT_DOT

    SMAT_DOT= matrix( (/ 0.0_dp,  0.0_dp,  -c2*dtheta(2), &
         & 0.0_dp, -s1*dtheta(1), c1*c2*dtheta(1)-s1*s2*dtheta(2),&
         & 0.0_dp, -c1*dtheta(1), -s1*c2*dtheta(1)-c1*s2*dtheta(2)/))

  end function find_angrate_dot_cosines


  ! ******************************************************************!
  ! Returns the inverse of the angular rate matrix for the supplied
  ! theta vector
  ! ******************************************************************!
  function find_angrate_inv_vec(thetain) &
       &result(SMAT_INV)

    type(vector) :: thetain
    type(matrix) :: SMAT_INV
    real(dp)     :: theta(NUM_SPAT_DIM)

    ! decompose the vector into array
    theta = array(thetain)

    ! call the method that takes array as input
    SMAT_INV =  find_angrate_inv_array(theta)

  end function find_angrate_inv_vec

  ! ******************************************************************!
  ! Returns the inverse of the angular rate matrix for the supplied
  ! theta array
  ! ******************************************************************!
  function find_angrate_inv_array(theta) &
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

    SMAT_INV =  find_angrate_inv_cosines(  c1, c2, c3, s1, s2, s3)

  end function find_angrate_inv_array

  ! ******************************************************************!
  ! Returns the inverse of the angular rate matrix for the supplied
  ! direction cosines
  ! ******************************************************************!
  function find_angrate_inv_cosines( c1, c2, c3, s1, s2, s3) &
       &result(SMAT_INV)

    real(dp)     :: c1, c2, c3, s1, s2, s3
    type(matrix) :: SMAT_INV

    SMAT_INV= matrix( (/  1.0_dp, s1*s2/c2, c1*s2/c2, &
         &0.0_dp, c1, -s1, 0.0_dp,&
         & s1/c2, c1/c2 /))

  end function find_angrate_inv_cosines

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
