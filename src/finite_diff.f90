!*********************************************************************!
! Module that performs finite differentiation
! Works for both scalar valued and vector valued functions
!---------------------------------------------------------------------!
! Scalar valued function:
!---------------------------------------------------------------------!
! R(q) 
!
! Example: 
!
! The function, R = sin(q) 
! The derivative: \frac{\partial R}{\partial q} = (R(q+h)-R(q))/h
!
! Here q is the variable
!---------------------------------------------------------------------!
! Vector valued function:
!---------------------------------------------------------------------!
! R(q, qdot, qdotdot ) 
!
! Example:
!
! R  = {R1; R2; R3}
! R1 = sin(q(1) *sin(q2) *sin(q(3) + c
! R2 = cos(q(1) *sin(q(3) + d 
! R3 = a*x+ b*y
! Here q(1), q(2), q(3) are the variables
! 
! The partial derivative of R with respect to q would be a 3x3 matrix
! In general the matrix size would be: [neq x nvar]
!
!********************************************************************!
module finite_diff

  implicit none

!!$  private
!!$
!!$  public
!!$
!!$  interface
!!$
!!$  end interface

  !-------------------------------------------------------------------!
  !  Define constants to manage precision [TUNABLE]
  !-------------------------------------------------------------------!
  integer, parameter      :: sp = kind(0.0)    ! single precision
  integer(sp), parameter  :: dp = kind(0.0d0)  ! double precision

  real(dp) :: h = 1.0d-10 ! step size for FD

  interface
     ! Signature of the function that the user must implement
!!$     function R(q)
!!$       integer, parameter      :: sp = kind(0.0)    ! single precision
!!$       integer(sp), parameter  :: dp = kind(0.0d0)  ! double precision
!!$       real(dp), dimension(:), intent(in) :: q
!!$       real(dp), dimension(size(q)) :: R
!!$     end function R
  end interface

contains

  !-------------------------------------------------------------------!
  ! The governing equation or function to be differentiated
  !-------------------------------------------------------------------!
  function FD(q, dh)

    integer, parameter                    :: sp = kind(0.0)    ! single precision
    integer(sp), parameter                :: dp = kind(0.0d0)  ! double precision

    real(dp), dimension(:)                :: q            ! variables 
    real(dp), dimension(size(q))          :: qtmp         ! temp storage

    
    real(dp), dimension(size(q), size(q)) :: FD           ! jacobian
    real(dp), optional                    :: dh           ! step size
    real(dp), dimension(size(q))          :: F, Ftmp, dF  ! solution vector   
    integer(sp)                           :: k, N

    ! R is to be implemented by the using program
    real(dp), external           :: R

    if (present(dh)) then
       write(*,*) "Using the supplied time step: ", dh
    else
       dh = 1.0e-10_dp
       write(*,*) "Using the default time step: ", dh
    end if

    ! compute original solution and perform init tasks
    F = R(q); qtmp = q; dF = 0.0_dp;

    N = size(q); 
    var: do  k = 1 , N ! loop through all variables     
       ! perturb k-th variable
       qtmp( k ) = qtmp (k)  + dh
       ! Find the perturbed solution
       Ftmp       = R(q)
       ! derivative of R wrt the k-th variables: partial R/partial q_k
       dF =  (Ftmp - F)/dh
       ! store in the jacobian matrix
       FD (:, k) = dF
       ! reset the variable to original state 
       ! and get ready to perturb the next variable
       qtmp = q
    end do var

  end function FD

end module finite_diff
