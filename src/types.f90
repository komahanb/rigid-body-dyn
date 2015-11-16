!*********************************************************************!
! Module that contains all the user defined data types
!
! Any new data types that are defined go here
!
!*********************************************************************!
module types

  use global_constants, only: sp, dp, NUM_SPAT_DIM

  implicit none

  !*******************************************************************!
  ! VECTOR datatype can be used for arrays that are represented in 
  ! three spatial dimensions, for example gravity, acceleration, 
  ! velocity, position, orientation etc
  !*******************************************************************!
  type vector
     real(dp)    :: x(NUM_SPAT_DIM) = 0.0_dp
  end type vector

  !*******************************************************************!
  ! MATRIX datatype can be used for matrices involved within spatial 
  ! dimensions, for example moment of inertia (J), rotation matrix(C, 
  ! Cdot), angular rates (S, Sdot)
  !*******************************************************************!
  type matrix
     real(dp)    :: ij(NUM_SPAT_DIM, NUM_SPAT_DIM) = 0.0_dp
  end type matrix

  ! ******************************************************************!
  ! JOINT datatype  characterizes the joints between two bodies
  ! ******************************************************************!
  type joint

     ! spherical, revolute, prismatic, planar
     character(len=10) :: type              

     ! the two interacting bodies
  !   type(body)        :: a , b             

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
