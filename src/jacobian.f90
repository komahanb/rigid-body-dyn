!=====================================================================!
! Module that assembles the jacobian
!---------------------------------------------------------------------!
! Terms in the jacobian come from 
!     (a) body 
!     (b) joints
!---------------------------------------------------------------------!
! Contains the function that assemble the below blocks into the full 
! jacobian
!---------------------------------------------------------------------!
!!$   type(matrix) :: D_R (4,4)
!!$   type(matrix) :: S_R (4,4)
!!$   type(matrix) :: S_S(1,1)
!!$   type(matrix) :: S_RS(4,1)
!!$   type(matrix) :: S_SR(1,4)
!=====================================================================!
module jacobian

  use global_variables, only: aa, bb
  use types, only: sp, dp, body, joint
  use utils

  implicit none

  ! all functions are private by default
  private

  ! only required functions are exposed to the calling program
  public get_jacobian
 
  !-------------------------------------------------------------------!
  ! The calling program will use get_jacobian as the interface to get
  ! the assembled jacobian
  !-------------------------------------------------------------------!
  ! Example usage:
  ! jac = get_jacobian(body) for a single body
  !  (only rigid motion terms are implemented)
  ! 
  ! jac = get_jacobian(body_array, joint_array) for multiple bodies 
  !   (not impl yet)
  !-------------------------------------------------------------------!
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

  !*******************************************************************!
  ! Jacobian of the equation of motion for the supplied BODY alpha
  !*******************************************************************!
  function D_R(alpha)

    type(body)    :: alpha
    type(matrix)  :: D_R(4,4) ! where 4 is the state vector length
    type(matrix)  :: O ! zero matrix
    type(matrix)  :: I ! identity matrix

    ! create local variables (for defn of these variables look 
    ! into types.f90)
    type(vector) :: r, theta, v, omega
    type(vector) :: r_dot,  theta_dot, v_dot, omega_dot
    type(vector) :: qs, qs_dot, qs_double_dot
    type(vector) :: c, fr, gr, f 
    type(matrix) :: J, p, h, K, M, C_mat, S, S_dot
    real(dp)     :: mass

    !--------------------------------------
    ! set the values for local variables
    !--------------------------------------

    ! rigid body state variables
    r             = alpha%r
    theta         = alpha%theta
    v             = alpha%v
    omega         = alpha%omega

    ! time derivative of states
    r_dot         = alpha%r_dot
    theta_dot     = alpha%theta_dot
    v_dot         = alpha%v_dot
    omega_dot     = alpha%omega_dot

    ! elatic state variables
    qs            = alpha%qs
    qs_dot        = alpha%qs_dot
    qs_double_dot = alpha%qs_double_dot

    ! other body Attributes
    mass          = alpha%mass
    c             = alpha%c
    J             = alpha%J

    p             = alpha%p
    h             = alpha%h

    K             = alpha%K
    M             = alpha%M

    C_mat         = alpha%C_mat
    S             = alpha%S
    S_dot         = alpha%S_dot

    fr            = alpha%fr
    gr            = alpha%gr

    f             = alpha%f

    ! some useful matrices
    O = matrix(zeros(NUM_SPAT_DIM))
    !    U = matrix(ones(NUM_SPAT_DIM))
    I = matrix(eye(NUM_SPAT_DIM))

    ! assemble jacobian below
    D_R(1,1) = aa * C_mat
    D_R(2,1) = O
    D_R(3,1) = O
    D_R(4,1) = O

    D_R(1,2) = skew(C_mat * r_dot) * S
    D_R(2,2) = S_dot + skew(S*theta_dot)*S + aa * S
    D_R(3,2) = O
    D_R(4,2) = O

    D_R(1,3) = -1.0_dp*I
    D_R(2,3) = O
    D_R(3,3) = mass*(aa*I + skew(omega))
    D_R(4,3) = aa*skew(c) + skew(c) * skew(omega)

    D_R(1,4) = O
    D_R(2,4) = -1.0_dp*I
    D_R(3,4) = -aa*skew(c) + skew(skew(c)*omega) - mass*skew(v) &
         &- skew(omega)*skew(c) 
    D_R(4,4) = aa*J  + skew(omega)*J - skew(J*omega) -skew(c)*skew(v)

  end function D_R

  !******************************************************************!
  ! Jacobian of the structural degree of freedom (S_R block)
  !*******************************************************************!
  function S_R(alpha)

    type(body)   :: alpha
    type(matrix) :: S_R(4,4) ! where 4 is the state vector length
    type(matrix)  :: O ! zero matrix

    O = matrix(zeros(NUM_SPAT_DIM))

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

  !*******************************************************************!
  ! Linear combination of stiffness and mass matrix S_S block
  !*******************************************************************!
  function S_S(alpha)

    type(body)   :: alpha
    type(matrix) :: S_S
    
    S_S = alpha%K  + bb*alpha%M ! K + b M

  end function S_S

  !*******************************************************************!
  ! Supplies the S_RS block in the jacobian for the supplies body
  !*******************************************************************!
  function S_RS(alpha)

    type(body)   :: alpha
    type(matrix) :: S_RS(1,4)
    type(matrix) :: O ! zero matrix

    O = matrix(zeros(NUM_SPAT_DIM))

    S_RS(1,1) = O
    S_RS(1,2) = O
    S_RS(1,3) = bb*alpha%p + aa*skew(alpha%omega)*alpha%p !bp +a wx p
    S_RS(1,4) = bb*alpha%h + aa*(skew(alpha%v)*alpha%p  &
         &+ skew(alpha%omega)*alpha%h) !bh +a(vx p + wx h)

  end function S_RS

  !*******************************************************************!
  ! Supplies the S_SR block in the jacobian for the supplies body
  !*******************************************************************!
  function S_SR(alpha)

    type(body)   :: alpha
    type(matrix) :: S_SR(4,1)
    type(matrix) :: O ! zero matrix

    O = matrix(zeros(NUM_SPAT_DIM))

    S_SR(1,1) = O
    S_SR(2,1) = O
    S_SR(3,1) = aa * trans(alpha%p) !a p^T
    S_SR(4,1) = aa * trans(alpha%p) !a h^T

  end function S_SR

end module jacobian
