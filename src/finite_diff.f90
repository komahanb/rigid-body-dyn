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
!---------------------------------------------------------------------!
! Note:
! 
!  The approximation of derivatives with finite-difference will cost
!  (N + 1) solutions of the governing equation which may not be 
!  readily available. 
! 
!  Also FD approximations are dependent on the step size chosen
!
!********************************************************************!
module finite_diff
  
  implicit none

  ! Use for systems like R(q) = 0
  interface finite_difference
     module procedure finite_diff_vector
  end interface finite_difference

  ! Use for systems like R(q,qdot) = 0
  interface finite_difference2
     module procedure finite_diff_vector2
  end interface finite_difference2

contains
  
  !*******************************************************************!
  ! The governing equation or function to be differentiated
  !-------------------------------------------------------------------!
  ! Input:
  ! q, qdot: the state variable and its time derivatives
  ! a : coefficient for linear sum of partial jacobians
  ! dh_in: the finite difference step size, optional
  !-------------------------------------------------------------------! 
  ! Output:
  ! jac: the jacobian computed with finite differemncing
  !*******************************************************************!
  function finite_diff_vector2(q, qdot, a, dh_in) result(jac)

    ! single precision
    integer, parameter                    :: sp = kind(0.0)   
    ! double precision
   integer(sp), parameter                 :: dp = kind(0.0d0)  
    ! input variables q, qdot
    real(dp), dimension(:)                :: q, qdot
    ! input step-size
    real(dp), optional                    :: dh_in      ! step size
    ! coefficient
    real(dp), intent(in) :: a 
    !returned jacobian
    real(dp), dimension(size(q), size(q)) :: jac
    ! temp jacobians
    real(dp), dimension(size(q), size(q)) :: jac1, jac2
    ! output
    real(dp), dimension(size(q))          :: F, Ftmp, dF, qtmp
    integer(sp)                           :: k, N
    real(dp)                              :: dh = 1.0e-10_dp !   size

    external                              :: residual2
     !-----------------------------------------------------------------!
    !--------------------- finite difference step --------------------!
    !-----------------------------------------------------------------!
    if (present(dh_in)) then ! if the step is supplied
       dh = dh_in
       !write(*,*) "Using the supplied step-size: ", dh
    else ! if the step is not supplied
       dh = 1.0e-10_dp
       !write(*,*) "Using the default time step-size: ", dh
    end if

    if (dh.eq.0.0_dp) stop " >> Wrong FD step-size! Please check!"

    !-----------------------------------------------------------------!
    !--------------- compute original sol (unperturbed)---------------!
    !-----------------------------------------------------------------!
    call residual2(q, qdot, F)

    !-----------------------------------------------------------------!
    !--------------- compute perturbed sol----------------------------!
    !-----------------------------------------------------------------!
    qtmp = q;     dF = 0.0_dp;     N = size(q); 
    state: do  k = 1 , N ! loop through all variables 
       ! N residual calls are made in this loop
       ! perturb k-th variable
       qtmp( k ) = qtmp (k)  + dh
       ! Find the perturbed solution by having q_dot fixed
       call residual2(qtmp, qdot, Ftmp)
       ! derivative of R wrt the k-th variables: partial R/partial q_k
       dF =  (Ftmp - F)/dh
       ! store in the jacobian matrix
       jac1 (:, k) = dF
       ! reset the variable to original state  and get ready to perturb the next variable
       qtmp = q
    end do state

    !-----------------------------------------------------------------!
    !--------------- compute perturbed sol----------------------------!
    !-----------------------------------------------------------------!
    qtmp = qdot;     dF = 0.0_dp;     N = size(qdot); 
    state_dot: do  k = 1 , N ! loop through all variables 
       ! N residual calls are made in this loop
       ! perturb k-th variable
       qtmp( k ) = qtmp (k)  + dh
       ! Find the perturbed solution
       call residual2(q, qtmp, Ftmp)
       ! derivative of R wrt the k-th variables: partial R/partial q_k
       dF =  (Ftmp - F)/dh
       ! store in the jacobian matrix
       jac2 (:, k) = dF
       ! reset the variable to original state  and get ready to perturb the next variable
       qtmp = qdot
    end do state_dot    

    ! linear sum of jacobians
    jac = jac1+a*jac2
    
  end function finite_diff_vector2

  !-------------------------------------------------------------------!
  ! The governing equation or function to be differentiated
  ! Input:
  ! The variables q as an array
  ! The finite difference step size, optional
  !-------------------------------------------------------------------!
  function finite_diff_vector(q, dh_in) result(jac)

    integer, parameter                    :: sp = kind(0.0)    ! single precision
    integer(sp), parameter                :: dp = kind(0.0d0)  ! double precision

    ! input variables
    real(dp), dimension(:)                :: q            ! variables 
    real(dp), optional                    :: dh_in           ! step size

    ! output
    real(dp), dimension(size(q), size(q)) :: jac           ! jacobian
    real(dp), dimension(size(q))          :: F, Ftmp, dF, qtmp
    integer(sp)                           :: k, N

    real(dp)                              :: dh = 1.0e-10_dp !   size

    ! R is to be implemented by the using program
    external                              :: residual

    !-----------------------------------------------------------------!
    !--------------------- finite difference step --------------------!
    !-----------------------------------------------------------------!
    if (present(dh_in)) then ! if the step is supplied
       dh = dh_in
       write(*,*) "Using the supplied step-size: ", dh
    else ! if the step is not supplied
       dh = 1.0e-10_dp
       write(*,*) "Using the default time step-size: ", dh
    end if

    if (dh.eq.0.0_dp) stop " >> Wrong FD step-size! Please check!"

    !-----------------------------------------------------------------!
    !--------------- compute original sol (unperturbed)---------------!
    !-----------------------------------------------------------------!
    call residual(q, F)

    !-----------------------------------------------------------------!
    !--------------- compute perturbed sol----------------------------!
    !-----------------------------------------------------------------!
    qtmp = q;     dF = 0.0_dp;     N = size(q); 
    variables: do  k = 1 , N ! loop through all variables 
       ! N residual calls are made in this loop
       ! perturb k-th variable
       qtmp( k ) = qtmp (k)  + dh
       ! Find the perturbed solution
       call residual(qtmp, Ftmp)
       ! derivative of R wrt the k-th variables: partial R/partial q_k
       dF =  (Ftmp - F)/dh
       ! store in the jacobian matrix
       jac (:, k) = dF
       ! reset the variable to original state  and get ready to perturb the next variable
       qtmp = q
    end do variables

  end function finite_diff_vector

end module finite_diff
