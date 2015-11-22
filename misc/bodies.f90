!=====================================================================!
! Module which handles BODIES.
!---------------------------------------------------------------------!
! Has functions to:
! (a) create_body based on supplied state, mass, J, c etc
! (b) update the state variables alone during time-integration
! (d) compute the rotation and angular rate matrices
! (c) 'toString' like function to print the body props (state+attrs)
!---------------------------------------------------------------------!
! Note: 
!---------------------------------------------------------------------!
! Use the interface names to access the functions instead of internal
! functions whose signature may change as per needs.
!=====================================================================!
module body_class

  use global_constants
  use global_variables
  use utils
  use dispmodule, only : disp

  implicit none

  ! Default make everything private
  private 

  ! Expose only the needed variables and functions
  public body, get_body
  public get_rotation
  public get_angrate, get_angrate_dot, get_angrate_inv
  public print_body

  !*******************************************************************!
  ! BODY datatype can be used to fully characterize the STATE and 
  ! ATTRIBUTES of a dynamic body .
  ! A body object contains virtually everything about the body 
  !*******************************************************************!
  type body

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
     ! elatic state variables and time derivatives
     !----------------------------------------------------------------!
     type(vector) :: qs
     type(vector) :: qs_dot
     type(vector) :: qs_double_dot

     !----------------------------------------------------------------!
     ! Body Attributes
     !----------------------------------------------------------------!

     real(dp)     :: mass           ! mass (denoted m in paper)   

     !     The format for c is: (in body frame)
     !     c = [ cx,  cy,  cz ]
     type(vector) :: c              ! first moment of inertia

     !  The format for J is: (in body frame)
     !  J = [ Jxx,  Jxy,  Jxz ] = [ J[0],  J[1],  J[2] ]
     !  . = [    ,  Jyy,  Jyz ] = [     ,  J[3],  J[4] ]
     !  . = [    ,     ,  Jzz ] = [     ,      ,  J[5] ]
     type(matrix) :: J              ! second moment of inertia

     type(matrix) :: p              ! 
     type(matrix) :: h              ! 

     type(matrix) :: K              ! stiffness matrix
     type(matrix) :: M              ! mass matrix

     type(matrix) :: C_mat          ! rotation matrix
     type(matrix) :: S
     type(matrix) :: S_dot          ! transformation matrix

     type(vector) :: fr             ! external/reaction force
     type(vector) :: gr             ! external/reaction torque

     type(vector) :: f              ! elastic force

     type(vector) :: g              ! gravity vector in local frame

     real(dp)     :: KE             ! kinetic energy of the body
     real(dp)     :: PE             ! potential energy of the body

  end type body


  !*******************************************************************!
  ! A common interface for creating bodies and updating existing ones
  !-------------------------------------------------------------------!
  ! The inputs are described in the methods under the interface. One 
  ! can use the same interface for updating existing bodies.
  ! Shortly: use this to create and update bodies
  !*******************************************************************!
  interface get_body
     module procedure get_rigid_body
  end interface get_body
  
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
  subroutine get_rigid_body(alpha, mass, c, J, fr, gr, q, qdot)

    ! inputs
    real(dp), optional, intent(in) :: mass 
    real(dp), optional, intent(in) :: c(NUM_SPAT_DIM)
    real(dp), optional, intent(in) :: J(NUM_SPAT_DIM, NUM_SPAT_DIM)
    real(dp), optional, intent(in) :: fr(NUM_SPAT_DIM)
    real(dp), optional, intent(in) :: gr(NUM_SPAT_DIM)

    real(dp), optional, intent(in) :: q(NDOF_PBODY)
    real(dp), optional, intent(in) :: qdot(NDOF_PBODY)

    ! input/output
    type(body), intent(inout) :: alpha

    !-----------------------------------------------------------------!
    ! Inertial properties of the body
    !-----------------------------------------------------------------!

    ! mass of the body
    if (present(mass)) alpha%mass = mass

    ! moment of inertia in body-fixed frame
    if (present(J)) alpha%J = matrix(J)       !-mass*skew(re)*skew(re)

    !first moment of inertia in body-fixed frame: mass*(cg location)
    if (present(c)) alpha%c = vector(c)       ! mass*re

    !-----------------------------------------------------------------!
    ! set the state into the body
    !-----------------------------------------------------------------!
    if (present(q)) then
       alpha%r         = vector(q(1:3))
       alpha%theta     = vector(q(4:6))
       alpha%v         = vector(q(7:9))
       alpha%omega     = vector(q(10:12))
    end if

    !-----------------------------------------------------------------!
    ! set the time derivatives of state into the body
    !-----------------------------------------------------------------!

    if (present(qdot)) then
       alpha%r_dot     = vector(qdot(1:3))
       alpha%theta_dot = vector(qdot(4:6))
       alpha%v_dot     = vector(qdot(7:9))
       alpha%omega_dot = vector(qdot(10:12))
    end if

    !-----------------------------------------------------------------!
    ! update the rotation and angular rate matrices
    !-----------------------------------------------------------------!

    alpha%C_mat     = get_rotation(alpha%theta)
    alpha%S         = get_angrate(alpha%theta)
    alpha%S_dot     = get_angrate_dot(alpha%theta, alpha%theta_dot)

    !-----------------------------------------------------------------!
    ! update the new direction of the gravity vector in body frame
    !-----------------------------------------------------------------!

    alpha%g = alpha%C_mat*GRAV

    !-----------------------------------------------------------------!
    ! Mechanical Energy 
    !-----------------------------------------------------------------!

    ! update the kinetic energy
    alpha%KE = 0.5_dp*(alpha%mass * alpha%v *  alpha%v &
         & + alpha%omega*alpha%J* alpha%omega) ! + coupling term

    ! update potential energy
    alpha%PE = alpha%mass*alpha%g*alpha%r

    !-----------------------------------------------------------------!
    ! Joint reactions
    !-----------------------------------------------------------------!

    ! reaction force 
    if (present(fr)) alpha%fr    = vector(fr)   

    ! reaction torque
    if (present(gr)) alpha%gr    = vector(gr)

    !call print_body(alpha)

  end subroutine get_rigid_body

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
  !* *****************************************************************!
  subroutine print_body(alpha)

    use dispmodule

    type(body):: alpha
    call disp('======================================================')
    call disp('---------------------BODY----------------------------')
    call disp('======================================================')
    call disp('')
    call disp('   > r          =   ', array(alpha%r), SEP=', ', &
         &ORIENT = 'ROW')
    call disp('   > theta      =   ', array(alpha%theta), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('   > v          =   ', array(alpha%v), SEP=', ', &
         &ORIENT = 'ROW')
    call disp('   > omega      =   ', array(alpha%omega), SEP=', ',&
         & ORIENT = 'ROW')
    call DISP('')
    call disp('   > r_dot      =   ', array(alpha%r_dot), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('   > theta_dot  =   ', array(alpha%theta_dot), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('   > v_dot      =   ', array(alpha%v_dot), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('   > omega_dot  =   ', array(alpha%omega_dot), SEP=', ',&
         &ORIENT = 'ROW')
    call DISP('')
    call disp('   > c          =   ', array(alpha%c), SEP=', ',&
         & ORIENT = 'ROW')
    call DISP('')
    call disp('   > J          =   ', matrix(alpha%J))
    call DISP('')
    call disp('   > Rot Mat    =   ', matrix(alpha%C_mat))
    call DISP('')
    call disp('   > Angrt Mat  =   ', matrix(alpha%S))
    call DISP('')
    call disp('   > S_dot      =   ', matrix(alpha%S_dot))
    call DISP('')
    call disp('   > fr         =   ', array(alpha%fr), SEP=', ', &
         &ORIENT = 'ROW')

    call disp('   > gr         =   ', array(alpha%gr), SEP=', ', &
         &ORIENT = 'ROW')
    call disp('')
    call disp('   > gravity    =   ', array(alpha%g), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('')
    call disp('   > Pot Energy =   ', alpha%PE)
    call disp('   > Kin Energy =   ', alpha%KE)
    call disp('   > Tot Energy =   ', alpha%KE + alpha%PE)
    call disp('======================================================')
  end subroutine print_body

end module body_class

!!$ ! ********************************************************
!!$  ! routine that returns the rotation matrix (euler angles)
!!$  ! based on the angles

! compiler threw ambiguity 
! error so I have commented out this impl for now and will explore
! later
!!$  ! ********************************************************
!!$  function get_rotation_from_state_vector(qr) result(CMAT)
!!$
!!$    real(dp), intent(in)       :: qr(NUM_STATES_PER_BODY)    
!!$    type(matrix), intent(out)  :: CMAT
!!$
!!$    CMAT = get_rotation_from_angles(qr(4:6)) ! angles are stored in indices 4:6
!!$
!!$  end function get_rotation_from_state_vector
!!$
!!$ ! ******************************************************************!
!!$  ! Returns the time derivative of rotation matrix (euler angles)
!!$  ! ******************************************************************!
!!$  function get_rotation_dot(theta)
!!$
!!$    real(dp) :: theta(NUM_SPAT_DIM)
!!$    type(matrix) :: get_rotation_dot
!!$
!!$    stop "dummy impl"

!!$ end function get_rotation_dot
!!$ !*******************************************************************!
!!$  ! create a body using the supplied state and other parameters
!!$  !-------------------------------------------------------------------!
!!$  ! q, qdot: state vector and time derivatives
!!$  ! mass   : mass of the body
!!$  ! re     : interested point on the body measure in body frame
!!$  !*******************************************************************!
!!$  function construct1(mass, re, q, qdot) result(alpha)
!!$
!!$    ! inputs
!!$    real(dp), intent(in)        :: mass 
!!$    type(vector), intent(in)    :: re 
!!$    real(dp), intent(in)        :: q(NDOF_PBODY)
!!$    real(dp), intent(in)        :: qdot(NDOF_PBODY)
!!$    ! output
!!$    type(body)                  :: alpha
!!$
!!$    ! local variables
!!$    real(dp), dimension(NUM_SPAT_DIM):: theta
!!$
!!$    ! make a local copy
!!$    theta = q(4:6)
!!$
!!$    ! rotation and rate matrices
!!$    alpha%C_mat     = get_rotation(theta)
!!$    alpha%S         = get_angrate(theta)
!!$    alpha%S_dot     = get_angrate_dot(theta, qdot(4:6))
!!$
!!$    !qdot = qdotin
!!$    !qdot(1:3) = array(trans(get_rotation(theta)) * vector(q(7:9)))
!!$    !qdot(4:6) = array(get_angrate_inv(theta) * vector(q(10:12)))
!!$
!!$    !  sets the state of the body and its time derivatives
!!$    call set_state(alpha, q, qdot, alpha)
!!$
!!$    ! set the gravity vector in body frame
!!$    alpha%g = alpha%C_mat*GRAV
!!$
!!$    ! mass of the body
!!$    alpha%mass  = mass
!!$
!!$    ! moment of inertia in body-fixed frame
!!$    alpha%J = -mass*skew(re)*skew(re)
!!$
!!$    !first moment of inertia in body-fixed frame: mass*(cg location)
!!$    alpha%c  = mass*re
!!$
!!$    ! reaction force 
!!$    alpha%fr = zeroV !mass*GRAV !? check
!!$!    alpha%fr =  vector( (/ cos(q(10)*time ), cos(q(11)*time), sin(q(12)*time ) /) )
!!$
!!$    ! reaction torque
!!$    alpha%gr = zeroV !alpha%J*alpha%omega_dot !? check
!!$
!!$    ! Translational + rotational KE = 0.5mV^2 + 0.5*wJw
!!$    ! (*) between two vectors does a dot product (see utils.f90)
!!$    ! Assuming the body axis to be at the centre of mass of the body
!!$    ! Inertial or body frame?
!!$    alpha%KE = 0.5_dp*(mass * alpha%v *  alpha%v &
!!$         & + alpha%omega*alpha%J* alpha%omega) ! + coupling term
!!$
!!$    ! Potential energy due to position. Should include strain energy later
!!$    ! Inertial or body frame?
!!$    alpha%PE = mass*alpha%g*alpha%r
!!$
!!$    !call print_body(alpha)
!!$
!!$  end function construct1


  !*******************************************************************!
  ! A common interface for setting (updating the state vectors) during
  ! time-integration. During time-integration only the state of the 
  ! bodies change and it is unnecessary to create new bodies for the 
  ! new state. Therefore the user can call this interface to update the 
  ! state of the supplied  bodies
  !-------------------------------------------------------------------! 
!!$  interface set_state
!!$     module procedure set_body_state
!!$  end interface set_state