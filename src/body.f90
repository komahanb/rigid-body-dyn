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
     
     real(dp)     :: m          ! mass

     type(vector) :: r, r_dot              ! radius of origin
     type(vector) :: theta, theta_dot      ! orientation of the body frame with respect to inertial

     type(vector) :: v, v_dot   ! the velocity and acceleration of the origin
     type(vector) :: omega, omega_dot    ! the angular velocity and acceleration of the body axis
     type(vector) :: qs, qs_dot, qs_double_dot  ! elastic state vectors due to deformation

     type(vector) :: F        ! external forces
     type(vector) :: G        ! external moments (torque)

     type(matrix) :: C        ! rotation matrix
     type(matrix) :: S        ! transformation matrix

     type(matrix) :: C_mat      ! first moment of inertia
     type(matrix) :: J        ! second moment of inertia
     type(matrix) :: K        ! stiffness matrix
     type(matrix) :: P        ! 
     real(dp)     :: a, b                       ! constant multipliers

  end type body
  

contains

! find the jacobian of the equation of motion for the supplied body
function jac(B)

  type(body)   :: B
  type(matrix) :: jac(4,4) ! where 4 is the state vector length

!  real(dp)     :: OO(num_spat_dim,num_spat_dim) ! zero matrix
!  real(dp)     :: UU(num_spat_dim,num_spat_dim) ! unit matrix
!  real(dp)     :: II(num_spat_dim,num_spat_dim) ! identity matrix
!  type(vector) :: vec


  type(matrix)  :: O(num_spat_dim,num_spat_dim) ! zero matrix
  type(matrix)  :: U(num_spat_dim,num_spat_dim) ! unit matrix
  type(matrix)  :: I(num_spat_dim,num_spat_dim) ! identity matrix

  O = matrix(zeros(num_spat_dim))
  U = matrix(ones(num_spat_dim))
  I = matrix(eye(num_spat_dim))
  
!!$  jac(1,1) = B%a * B%C_mat
!!$  jac(2,1) = skew(B%C_mat * B%r_dot) * B%S
!!$  jac(3,1) = -1._dp*U
!!$  jac(4,1) = O
!!$
!!$  jac(1,2) = O
!!$  jac(2,2) = B%S_dot + skew(B%S*B%theta_dot)*B%S + B%a * B%S
!!$  jac(3,2) = O
!!$  jac(4,2) = -1.0_dp*UU
!!$
!!$  jac(1,3) = O
!!$  jac(2,3) = O
!!$  jac(3,3) = B%m*(B%a*U + skew(B%omega))
!!$  jac(4,3) = -B%a*skew(B%C) + skew(skew(B%C)*B%omega) - B%m*skew(B%V) - skew(B%omega)*skew(B%C)
!!$
!!$  jac(1,4) = O
!!$  jac(2,4) = O
!!$  jac(3,4) = B%a*skew(B%C) + skew(B%C) * skew(B%omega)
!!$  jac(4,4) = B%a*B%J  + skew(B%omega)*B%J - skew(B%J*B%omega) -skew(B%c)*skew(B%V)
!!$

!  vec = matrix_vector(B%C_mat,B%r_dot)
!  print*, skew() * B%S)

!  print*,matrix(skew(vec))
!  print*, matrix(B%S)
!  print*, matrix_matrix(skew(vec) ,B%S )
! print*, skew(vec) * B%S
!  print*, skew(B%C_mat * B%r_dot) * B%S !! correct as of now
!  print*, cross(B%C_mat * B%r_dot, B%S)
!  print*, B%C_mat * B%r_dot * B%S ! vector right now should be matrix

end function jac

end module test_body
