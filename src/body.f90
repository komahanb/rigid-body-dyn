module body_mod

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
     type(matrix) :: M_mat                      ! mass matrix

     type(matrix) :: p                          ! 
     type(matrix) :: h                          ! 

     real(dp)     :: a, b                       ! constant multipliers
     

     ! for elastic
     ! type(vector) :: q_dot !? maybe qs

  end type body
  
end module body_mod
