!*********************************************************************************
! Module that contains all the user defined data types that are used in the code
!
! Any new data types that are defined go here
!
!*********************************************************************************
module types

  use global_constants, only: sp, dp, NUM_SPAT_DIM

  implicit none

  !*********************************************************************************
  ! VECTOR datatype can be used for arrays that are represented in three spatial dimensions
  ! example gravity, acceleration, velocity, position, orientation etc
  !*********************************************************************************
  type vector
     real(dp)    :: x(NUM_SPAT_DIM) = 0.0_dp
  end type vector

  !*********************************************************************************
  ! MATRIX datatype can be used for matrices involved within spatial dimensions
  ! example moment of inertia (J), rotation matrix(C, Cdot), angular rates (S, Sdot)
  !*********************************************************************************
  type matrix
     real(dp)    :: ij(NUM_SPAT_DIM, NUM_SPAT_DIM) = 0.0_dp
  end type matrix

  !*********************************************************************************
  ! BODY datatype can be used to fully characterize the STATE and ATTRIBUTES of a dynamic body
  ! A body object contains virtually everything about the body 
  !*********************************************************************************
  type body

     ! rigid body state variables
     type(vector) :: r              ! origin of the body frame
     type(vector) :: theta          ! orientation of the body frame with inertial (euler angles)
     type(vector) :: v              ! velocity of the origin with respect to inertial
     type(vector) :: omega          ! angular velocity of the body frame with respect to inertial

     ! time derivative of states
     type(vector) :: r_dot
     type(vector) :: theta_dot
     type(vector) :: v_dot
     type(vector) :: omega_dot

     ! elatic state variables
     type(vector) :: qs
     type(vector) :: qs_dot
     type(vector) :: qs_double_dot

     ! other body Attributes
     real(dp)     :: mass           ! mass (denoted m in paper)   
     type(vector) :: c              ! first moment of inertia
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

  end type body

  ! ********************************************************************
  ! JOINT datatype is used to characterize the joints between two bodies
  ! ********************************************************************
  type joint

     character(len=10) :: type              ! spherical, revolute, prismatic, planar
     type(body)        :: a , b             ! the two interacting bodies
     type(vector)      :: aPoint, bPoint    ! point of interaction as measured in their body axes

  end type joint

end module types
