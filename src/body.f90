module test_body
  use utils
  use matrix_utils
  implicit none
  ! each body shall contain the following
  ! body fixed frame located at O position(R, Theta) wrt inertial
  ! the velocity v0 (body)
  ! the angular omega (body)
  ! A rotation matrix (inertia-body)
  ! A transformation matrix (inertial-body)

  ! mass moment of inertia
  ! stiffness matrix
  ! EQ1, EQ2, EQ3 = m{v0_dot} - {c} x {omega_dot}  + {p}.{qs_double_dot} + {omega} x m {v0} - {omega} x {c} x {omega} + {omega} x {p}.{qs_dot}
  ! EQ4, EQ5, EQ6 = {c}.{v0_dot}  + [J]{omega_dot} + [h]{qs_double_dot} + [c] x {omega} x {v0} + {omega} x [J] {omega} + {}x[p].{qs_dot} + {omega} x [h]{qs_dot} =gr

  type body
     
     real(dp)     :: m                          ! mass

     type(vector) :: r, r_dot                   ! radius of origin
     type(vector) :: theta, theta_dot           ! orientation of the body frame with respect to inertial

     type(vector) :: v                          ! the velocity and acceleration of the origin
     type(vector) :: omega                      ! the angular velocity and acceleration of the body axis
     type(vector) :: qs, qs_dot, qs_double_dot  ! elastic state vectors due to deformation

     type(vector) :: F                          ! external forces
     type(vector) :: G                          ! external moments (torque)

     type(matrix) :: C_mat                      ! rotation matrix
     type(matrix) :: S, S_dot                   ! transformation matrix

     type(vector) :: c                          ! first moment of inertia
     type(matrix) :: J                          ! second moment of inertia
     type(matrix) :: K                          ! stiffness matrix

     type(matrix) :: p                          ! 
     type(matrix) :: h                          ! 

     real(dp)     :: a, b                       ! constant multipliers


     ! for elastic
     ! type(vector) :: q_dot !? maybe qs

  end type body
  

contains

! find the jacobian of the equation of motion for the supplied BODY alpha
function D_R(alpha)

  type(body)   :: alpha
  type(matrix) :: D_R(4,4) ! where 4 is the state vector length
  type(matrix)  :: O ! zero matrix
  type(matrix)  :: U ! unit matrix
  type(matrix)  :: I ! identity matrix

  ! some useful matrices
  O = matrix(zeros(num_spat_dim))
  U = matrix(ones(num_spat_dim))
  I = matrix(eye(num_spat_dim))
  
  ! assemble jacobian
  D_R(1,1) = alpha%a * alpha%C_mat
  D_R(2,1) = skew(alpha%C_mat * alpha%r_dot) * alpha%S
  D_R(3,1) = -1.0_dp*U
  D_R(4,1) = O

  D_R(1,2) = O
  D_R(2,2) = alpha%S_dot + skew(alpha%S*alpha%theta_dot)*alpha%S + alpha%a * alpha%S
  D_R(3,2) = O
  D_R(4,2) = -1.0_dp*U

  D_R(1,3) = O
  D_R(2,3) = O
  D_R(3,3) = alpha%m*(alpha%a*U + skew(alpha%omega))
  D_R(4,3) = -alpha%a*skew(alpha%c) + skew(skew(alpha%c)*alpha%omega)  - alpha%m*skew(alpha%V) - skew(alpha%omega)*skew(alpha%C)

  D_R(1,4) = O
  D_R(2,4) = O
  D_R(3,4) = alpha%a*skew(alpha%C) + skew(alpha%C) * skew(alpha%omega)
  D_R(4,4) = alpha%a*alpha%J  + skew(alpha%omega)*alpha%J - skew(alpha%J*alpha%omega) -skew(alpha%c)*skew(alpha%V)

!  real(dp)     :: OO(num_spat_dim,num_spat_dim) ! zero matrix
!  real(dp)     :: UU(num_spat_dim,num_spat_dim) ! unit matrix
!  real(dp)     :: II(num_spat_dim,num_spat_dim) ! identity matrix
!  type(vector) :: vec

!  vec = matrix_vector(B%C_mat,B%r_dot)
!  print*, skew() * B%S)

!  print*,matrix(skew(vec))
!  print*, matrix(B%S)
!  print*, matrix_matrix(skew(vec) ,B%S )
! print*, skew(vec) * B%S
!  print*, skew(B%C_mat * B%r_dot) * B%S !! correct as of now
!  print*, cross(B%C_mat * B%r_dot, B%S)
!  print*, B%C_mat * B%r_dot * B%S ! vector right now should be matrix

end function D_R


! jacobian of the structural degree of freedom
function S_R(alpha)

  type(body)   :: alpha
  type(matrix) :: S_R(4,4) ! where 4 is the state vector length
  type(matrix)  :: O ! zero matrix

  O = matrix(zeros(num_spat_dim))

  S_R(1,1) = O
  S_R(2,1) = O
  S_R(3,1) = O
  S_R(4,1) = O

  S_R(1,2) = O
  S_R(2,2) = O
  S_R(3,2) = O
  S_R(4,2) = O

  S_R(1,3) = O
  S_R(2,3) = O
  S_R(3,3) = O
  S_R(4,3) = - skew(alpha%p*alpha%qs_dot)

  S_R(1,4) = O
  S_R(2,4) = O
  S_R(3,4) = - skew(alpha%p*alpha%qs_dot)
  S_R(4,4) = - skew(alpha%h*alpha%qs_dot)

end function S_R


end module test_body
