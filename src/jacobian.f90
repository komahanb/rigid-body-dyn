
!=========================================================
! module that assembles the jacobian
!=========================================================
! Terms in the jacobian come from 
!     (a) body 
!     (b) joints
!=========================================================
! Contains the function that assemble the below matrices:
!=========================================================
!!$   type(matrix) :: D_R (4,4)
!!$   type(matrix) :: S_R (4,4)
!!$   type(matrix) :: S_S(1,1)
!!$   type(matrix) :: S_RS(4,1)
!!$   type(matrix) :: S_SR(1,4)
!=========================================================
module jacobian

  use global_variables, only: aa, bb
  use types, only: sp, dp, body, joint
  use utils

  implicit none

  ! all functions are private by default
  private
  
  ! only required functions are exposed to the calling program
  public get_jacobian
  
  ! The calling program will use get_jacobian as the interface to get the assembled jacobian
  interface get_jacobian
     module procedure jac_rigid
  end interface get_jacobian

contains

  !****************************************************
  ! returns the rigid body terms in the jacobian matrix
  !****************************************************
  function jac_rigid(alpha)

    type(body)      :: alpha
    real(dp)        :: jac_rigid(12,12)
    type(matrix)    :: DR1 (4,4)

    ! assemble the jacobian of type(matrix)
    DR1              = D_R(alpha)   

    ! now return in real-matrix form needed for computations
    jac_rigid       = matrix(DR1)

  end function jac_rigid

  !****************************************************
  ! returns the FULL jacobian matrix terms
  !****************************************************
  function jac_full(alpha) !?! may be more inputs e.g. joint

    type(body)      :: alpha
    real(dp)        :: jac_full(1,1) !?? decide the dim

    ! constituent sub-matrices in the jacobian
    type(matrix) :: DR (4,4)
    type(matrix) :: SR (4,4)
    type(matrix) :: SS(1,1)
    type(matrix) :: SRS(4,1)
    type(matrix) :: SSR(1,4)

    stop"not implemted yet"

  end function jac_full

  !*************************************************************************
  ! find the jacobian of the equation of motion for the supplied BODY alpha
  !*************************************************************************
  function D_R(alpha)
    use dispmodule
    type(body)    :: alpha
    type(matrix)  :: D_R(4,4) ! where 4 is the state vector length
    type(matrix)  :: O ! zero matrix
    !    type(matrix)  :: U ! unit matrix
    type(matrix)  :: I ! identity matrix

    ! some useful matrices
    O = matrix(zeros(num_spat_dim))
    !    U = matrix(ones(num_spat_dim))
    I = matrix(eye(num_spat_dim))

    ! assemble jacobian
    D_R(1,1) = aa * alpha%C_mat
    D_R(2,1) = O
    D_R(3,1) = O
    D_R(4,1) = O

    D_R(1,2) = skew(alpha%C_mat * alpha%r_dot) * alpha%S
    D_R(2,2) = alpha%S_dot + skew(alpha%S*alpha%theta_dot)*alpha%S + aa * alpha%S
    D_R(3,2) = O
    D_R(4,2) = O

    D_R(1,3) = -1.0_dp*I
    D_R(2,3) = O
    D_R(3,3) = alpha%mass*(aa*I + skew(alpha%omega))
    D_R(4,3) = aa*skew(alpha%c) + (skew(alpha%c) * skew(alpha%omega))

    D_R(1,4) = O
    D_R(2,4) = -1.0_dp*I
    D_R(3,4) = -aa*skew(alpha%c) + skew(skew(alpha%c)*alpha%omega)  - alpha%mass*skew(alpha%v) - skew(alpha%omega)*skew(alpha%c) 
    D_R(4,4) = aa*alpha%J  + skew(alpha%omega)*alpha%J - skew(alpha%J*alpha%omega) -skew(alpha%c)*skew(alpha%V)

  end function D_R

  !**********************************************
  ! jacobian of the structural degree of freedom
  !**********************************************
  function S_R(alpha)

    type(body)   :: alpha
    type(matrix) :: S_R(4,4) ! where 4 is the state vector length
    type(matrix)  :: O ! zero matrix

    O = matrix(zeros(num_spat_dim))

    S_R(1,1) = O
    S_R(1,2) = O
    S_R(1,3) = O
    S_R(1,4) = O

    S_R(2,1) = O
    S_R(2,2) = O
    S_R(2,3) = O
    S_R(2,4) = O

    S_R(3,1) = O
    S_R(3,2) = O
    S_R(3,3) = O
    S_R(3,4) = - skew(alpha%p*alpha%qs_dot)

    S_R(4,1) = O
    S_R(4,2) = O
    S_R(4,3) = - skew(alpha%p*alpha%qs_dot)
    S_R(4,4) = - skew(alpha%h*alpha%qs_dot)

  end function S_R

  !****************************************************
  ! linear combination of stiffness and mass matrix
  !****************************************************
  function S_S(alpha)

    type(body)   :: alpha
    type(matrix) :: S_S

    S_S = alpha%K  + bb*alpha%M ! K + b M

  end function S_S

  !****************************************************
  !****************************************************
  function S_RS(alpha)

    type(body)   :: alpha
    type(matrix) :: S_RS(1,4)
    type(matrix) :: O ! zero matrix

    O = matrix(zeros(num_spat_dim))

    S_RS(1,1) = O
    S_RS(1,2) = O
    S_RS(1,3) = bb*alpha%p + aa*skew(alpha%omega)*alpha%p !bp +a wx p
    S_RS(1,4) = bb*alpha%h + aa*(skew(alpha%v)*alpha%p  + skew(alpha%omega)*alpha%h) !bh +a(vx p + wx h)

  end function S_RS

  !****************************************************
  !****************************************************
  function S_SR(alpha)

    type(body)   :: alpha
    type(matrix) :: S_SR(4,1)
    type(matrix) :: O ! zero matrix

    O = matrix(zeros(num_spat_dim))

    S_SR(1,1) = O
    S_SR(2,1) = O
    S_SR(3,1) = aa * trans(alpha%p) !a p^T
    S_SR(4,1) = aa * trans(alpha%p) !a h^T

  end function S_SR

end module jacobian
