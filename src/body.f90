! module that defines the configuration of a single body
module body_mod
  use utils
  use matrix_utils
  implicit none

  private

  public body, CBI, SIB, alpha, R_rigid, setstatevars, create_body, print_body

  type body

     ! rigid body state variables
     type(vector) :: r
     type(vector) :: theta
     type(vector) :: v
     type(vector) :: omega

     type(vector) :: r_dot
     type(vector) :: theta_dot
     type(vector) :: v_dot
     type(vector) :: omega_dot

     ! elatic state variables
     type(vector) :: qs
     type(vector) :: qs_dot
     type(vector) :: qs_double_dot

     ! other body properties
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

  type(body)      :: alpha

contains

  subroutine print_body(alpha)
    use dispmodule
    implicit none

    type(body):: alpha

    CALL DISP('.....................Body.........................')
    call DISP('')
    CALL DISP('   r          =   ', get_array(alpha%r), SEP=', ', ORIENT = 'ROW')
    CALL DISP('   theta      =   ', get_array(alpha%theta), SEP=', ', ORIENT = 'ROW')
    CALL DISP('   v          =   ', get_array(alpha%r), SEP=', ', ORIENT = 'ROW')
    CALL DISP('   omega      =   ', get_array(alpha%theta), SEP=', ', ORIENT = 'ROW')
    call DISP('')
    CALL DISP('   r_dot      =   ', get_array(alpha%r), SEP=', ', ORIENT = 'ROW')
    CALL DISP('   theta_dot  =   ', get_array(alpha%theta), SEP=', ', ORIENT = 'ROW')
    CALL DISP('   v_dot      =   ', get_array(alpha%r), SEP=', ', ORIENT = 'ROW')
    CALL DISP('   omega_dot  =   ', get_array(alpha%theta), SEP=', ', ORIENT = 'ROW')
    call DISP('')
    CALL DISP('   c          =   ', get_array(alpha%c), SEP=', ', ORIENT = 'ROW')
    call DISP('')
    CALL DISP('   J          =   ', get_matrix(alpha%J))
    call DISP('')
    CALL DISP('   C          =   ', get_matrix(alpha%C_mat))
    call DISP('')
    CALL DISP('   S          =   ', get_matrix(alpha%S))
    call DISP('')
    CALL DISP('   S_dot      =   ', get_matrix(alpha%S_dot))
    call DISP('')
    CALL DISP('   fr         =   ', get_array(alpha%fr), SEP=', ', ORIENT = 'ROW')
    CALL DISP('   gr         =   ', get_array(alpha%gr), SEP=', ', ORIENT = 'ROW')
    call DISP('')
    CALL DISP('..................................................')

  end subroutine print_body

  ! assembles the R_rigid vector for the body
  function R_rigid(alpha)

    type(vector) :: R_rigid(4)
    type(body) :: alpha

    R_rigid(1)  = alpha%C_mat*alpha%r_dot - alpha%v

    R_rigid(2)  = alpha%S*alpha%theta_dot - alpha%omega

    R_rigid(3)  = alpha%mass*alpha%v_dot - skew(alpha%c)*alpha%omega_dot +alpha%p*alpha%qs_double_dot &
         &+ skew(alpha%omega)*(alpha%mass*alpha%v - alpha%c*alpha%omega + alpha%p*alpha%qs_dot) - alpha%fr

    R_rigid(4)  = skew(alpha%c)*alpha%v_dot + alpha%J*alpha%omega_dot + alpha%h*alpha%qs_double_dot &
         &+ skew(alpha%c)*skew(alpha%omega)*alpha%v + skew(alpha%omega)*alpha%J*alpha%omega &
         &+ skew(alpha%v)*alpha%p*alpha%qs_dot + skew(alpha%omega)*alpha%h*alpha%qs_dot &
         &+ skew(alpha%omega)*alpha%h*alpha%qs_dot -alpha%gr

  end function R_rigid


  ! assembles the R_elastic vector for the body
  function R_elastic(alpha)
    type(vector) :: R_elastic
    type(body) :: alpha
    stop"dummy impl"
  end function R_elastic


  ! function to create body and set the class variables 
!!$function makeBody()
!!$type(body) :: makebody
!!$end function makeBody



  ! sets the 
  subroutine SetStateVars(y, yprime, alpha)
    !,  r, theta, v, omega, r_dot, theta_dot, v_dot, omega_dot)
    !  use constants
    !  use utils
    !  use body_mod,only:alpha
    implicit none

    real(dp), intent(in)          :: y(12), yprime(12)
    real(dp)                      :: r(num_spat_dim), theta(num_spat_dim), v(num_spat_dim), omega(num_spat_dim)
    real(dp)                      :: r_dot(num_spat_dim), theta_dot(num_spat_dim), v_dot(num_spat_dim), omega_dot(num_spat_dim)
    type(body),intent(inout)      :: alpha

!!$  r       = (/ 0.0_dp, 0._dp, 0.0_dp /)
!!$  theta   = (/ deg2rad(10.0d0), deg2rad(0.0d0), deg2rad(0.0d0) /)
!!$  v     = (/ 1.0d0, 0.0d0, 0.0d0 /)
!!$  omega = (/ 1.0d0, 0.0d0, 0.0d0 /)

    alpha%r         = vector(Y(1:3))
    alpha%theta     = vector(Y(4:6))
    alpha%v         = vector(Y(7:9))
    alpha%omega     = vector(Y(10:12))

    ! define velocities of the body frame (q_dot terms)
    alpha%r_dot     = vector(YPRIME(1:3))
    alpha%theta_dot = vector(YPRIME(4:6))
    alpha%v_dot     = vector(YPRIME(7:9))
    alpha%omega_dot = vector(YPRIME(10:12))

  end subroutine SetStateVars

  !**************************************************************
  ! create a body using the supplied state and other parameters
  !**************************************************************
  subroutine create_body(m, g0, re, Y, YPRIME, alpha)

    use constants
    use utils
    use matrix_utils

    implicit none

    real(dp), intent(in)        :: m 
    type(vector), intent(in)    :: re , g0
    real(dp), intent(in)        :: y(12), yprime(12)
    type(body), intent(out)     :: alpha
    real(dp)                    :: theta(3)

    call SetStateVars(Y, YPRIME, alpha)

    theta = Y(4:6)

    ! rotation and rate matrices
    alpha%C_mat     = matrix(CBI(theta)) !?
    alpha%S         = matrix(SIB(theta)) !?

    ! mass of the body
    alpha%mass  = m

    ! mass moment of intertia (second moment)
    !  temp     = (2.0_dp*m*abs(re)**2)/5.0_dp  
    !  alpha%J  = temp*idM
    alpha%J = -m*skew(re)*skew(re)

    ! first moment of mass
    alpha%c  = m*re

    ! force 
    alpha%fr = m*g0

    ! torque
    alpha%gr = alpha%J*alpha%omega_dot

    !  print*, "new body=", alpha

    call print_body(alpha)

  end subroutine create_body

end module body_mod

!!$
!!$!(r, theta, v, omega, r_dot, theta_dot, v_dot, omega_dot,&
!!$!     & mass, c, J, p, h, K, M, C_mat, S, S_dot, fr, gr)
!!$

! define positions of body frame from inertial (q terms)
!!$  r       = (/ 0.0_dp, 0._dp, 0.0_dp /)
!!$  theta   = (/ deg2rad(10.0d0), deg2rad(0.0d0), deg2rad(0.0d0) /)
!!$  v     = (/ 1.0d0, 0.0d0, 0.0d0 /)
!!$  omega = (/ 1.0d0, 0.0d0, 0.0d0 /)


!!$
!!$  alpha%r         = vector(r)
!!$  alpha%theta     = vector(theta)
!!$  alpha%v         = vector(v)
!!$  alpha%omega     = vector(omega)
!!$
!!$  ! define velocities of the body frame (q_dot terms)
!!$  alpha%r_dot     = zeroV
!!$  alpha%theta_dot = zeroV
!!$  alpha%v_dot     = zeroV
!!$  alpha%omega_dot = zeroV
