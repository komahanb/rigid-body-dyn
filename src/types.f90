!*********************************************************************!
! Module that contains all the user defined data types
!
! Any new data types that are defined go here
!
!*********************************************************************!
module types

  use global_constants, only: sp, dp, NUM_SPAT_DIM

  implicit none
  
  !-------------------------------------------------------------------!
  ! Type that is used for conversion between real and complex types
  !-------------------------------------------------------------------!
  type scalar
      real(dp)    :: x = 0.0_dp
!     real(dp)    :: y = 0.0_dp
!     complex(dp) :: z
  end type scalar
  
  !-------------------------------------------------------------------!
  ! VECTOR datatype can be used for arrays that are represented in 
  ! three spatial dimensions, for example gravity, acceleration, 
  ! velocity, position, orientation etc
  !-------------------------------------------------------------------!
  type vector
     real(dp)    :: x(NUM_SPAT_DIM) = 0.0_dp
  end type vector

  !-------------------------------------------------------------------!
  ! MATRIX datatype can be used for matrices involved within spatial 
  ! dimensions, for example moment of inertia (J), rotation matrix(C, 
  ! Cdot), angular rates (S, Sdot)
  !-------------------------------------------------------------------!
  type matrix
     real(dp)    :: ij(NUM_SPAT_DIM, NUM_SPAT_DIM) = 0.0_dp
  end type matrix

  !-------------------------------------------------------------------!
  ! BODY datatype can be used to fully characterize the STATE and 
  ! ATTRIBUTES of a dynamic body .
  ! A body object contains virtually everything about the body 
  !-------------------------------------------------------------------!
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

  ! ------------------------------------------------------------------!
  ! JOINT datatype  characterizes the joints between two bodies
  ! ------------------------------------------------------------------!
  type joint

     ! spherical, revolute, prismatic, planar
     character(len=10) :: type              

     ! the two interacting bodies
     type(body)        :: a , b             

     ! point of interaction as measured in their body axes
     type(vector)      :: aPoint, bPoint    

  end type joint


  !-------------------------------------------------------------------!
  ! Type that stores the jacobian in sparse matrix format
  !-------------------------------------------------------------------!
  type jac_mat
     integer(sp) :: n
     integer(sp) :: nz
     real(dp),allocatable :: row(:)
     real(dp),allocatable :: col(:)
     real(dp),allocatable :: val(:)
  end type jac_mat

end module types
