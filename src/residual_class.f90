!=====================================================================!
! Module that assembles the residual
!---------------------------------------------------------------------!
! The functions to call are by the name of the interface. 
!
! It is recommended to refrain from calling by name of the internal 
! functions directly as their signatures may change in future to adapt
! different requirements and scenarios.
!
!---------------------------------------------------------------------!
! This module is implemented with the following preference in mind.
!
! Prefers passing derived datatypes such as BODY and JOINT instead of 
! data in primitive form such as the arrays (Q, QDOT) because:
!
! (a) it is cleaner to pass data in and out as BEANS compared to
! primitive datatypes --that are more susceptible to programmer errors
! like passing the wrong arguments
!
! (b) future enhancements come easy when data is passed in BEAN form 
! and can be easily done without having to change the function 
! signatures. One can just add extra variables into the BODY or JOINT
! types instead of sending one or more variables as new arguments.

! For example: When extending from rigid-case to elastic case, the 
! second derivatives of the state (q_double_dot) are needed to assemble
! the residual, which may warrant a function signature change when 
! data is passed in an out in primitive form.
! 
! (c) For multi-body systems it is cleaner to keep track of bodies and 
! joints, compared to a very long Q, Q_dot, Q_double_dot arrays whose
! indices may correspond to any body or joint eqn.
!
!---------------------------------------------------------------------!
!  Example usage:
!  res(1:12) = get_residual(body)
!  res(:)    = get_residual(body_array, joint_array) (not fully func)
!---------------------------------------------------------------------!
! Inputs:
! array of bodies and joints (multi-body dynamics) OR (not fully func)
! single body (single body dynamics)
!---------------------------------------------------------------------!
! Output:
! Residual array of size NUM_STATES
!=====================================================================!
module residual_class

  use global_constants, only: dp
!!$  use global_variables
!!$  use utils, only: operator(*), operator(+), operator(-),&
!!$       & matrix, array, skew
!!$  use types, only: matrix, vector
!!$  use body_class, only : body
!!$  use rigid_body_class, only : rigid_body
!!$  use joint_class, only: joint

  implicit none

  ! Makes all functions private by default
  private

  ! Expose only the needed functions. 
  public residual
  
  type, abstract :: residual

     private

     real(dp), dimension(:), allocatable :: R

   contains

     procedure :: get_residual
     procedure :: set_residual

  end type residual

contains

  !*******************************************************************!
  ! Getter for the residual
  !*******************************************************************!
  
  function get_residual(this)

    class(residual) :: this
    real(dp), dimension(:), allocatable  :: get_residual
    
    allocate( get_residual( size(this % R) ) )

    get_residual = this % R
    
  end function get_residual

  !*******************************************************************!
  ! Setter for the residual
  !*******************************************************************!
  
  subroutine set_residual(this, R)

    class(residual) :: this
    real(dp), dimension(:)  :: R
    
    allocate(this % R (size (R)))

    this % R  =  R

  end subroutine set_residual

end module residual_class
