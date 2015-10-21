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

  private
  
  public body, CBI, SIB
  
  type body
     
     real(dp)     :: m                          ! mass

     type(vector) :: r
     type(vector) :: theta
     type(vector) :: v
     type(vector) :: omega

     type(vector) :: r_dot                   ! radius of origin
     type(vector) :: theta_dot           ! orientation of the body frame with respect to inertial
     type(vector) :: v_dot                          ! the velocity and acceleration of the origin
     type(vector) :: omega_dot                      ! the angular velocity and acceleration of the body axis

     type(vector) :: qs
     type(vector) :: qs_dot
     type(vector) :: qs_double_dot  ! elastic state vectors due to deformation

     type(vector) :: c                          ! first moment of inertia
     type(matrix) :: J                          ! second moment of inertia

     type(matrix) :: K                          ! stiffness matrix
     type(matrix) :: M_mat                      ! mass matrix

     type(matrix) :: p                          ! 
     type(matrix) :: h                          ! 

 
     type(matrix) :: C_mat                      ! rotation matrix
     type(matrix) :: S, S_dot                   ! transformation matrix

     type(vector) :: fr                          ! external forces
     type(vector) :: gr                          ! external moments (torque)
     
     ! for elastic
     ! type(vector) :: q_dot !? maybe qs

  end type body

contains

! residual vector
function residual(alpha)

  type(vector) :: residual(4)
  type(body) :: alpha
  
  residual(1)  = alpha%C_mat*alpha%r_dot - alpha%v

  residual(2)  = alpha%S*alpha%theta_dot - alpha%omega

  residual(3)  = alpha%m*alpha%v_dot - skew(alpha%c)*alpha%omega_dot +alpha%p*alpha%qs_double_dot &
       &+ skew(alpha%omega)*(alpha%m*alpha%v - alpha%c*alpha%omega + alpha%p*alpha%qs_dot) - alpha%fr

  residual(4)  = skew(alpha%c)*alpha%v_dot + alpha%J*alpha%omega_dot + alpha%h*alpha%qs_double_dot &
       &+ skew(alpha%c)*skew(alpha%omega)*alpha%v + skew(alpha%omega)*alpha%J*alpha%omega &
       &+ skew(alpha%v)*alpha%p*alpha%qs_dot + skew(alpha%omega)*alpha%h*alpha%qs_dot &
       &+ skew(alpha%omega)*alpha%h*alpha%qs_dot -alpha%gr

end function residual

! rotates a vector from body frame inertial
function CBI(theta)
  real(dp)                :: theta(num_spat_dim)
  real(dp)                :: CBI(num_spat_dim, num_spat_dim)

  CBI(1,1) =  cos(theta(2))*cos(theta(3)) + sin(theta(1))*sin(theta(2))*sin(theta(3))
  CBI(2,1) =  cos(theta(1))*sin(theta(3))
  CBI(3,1) = -sin(theta(2))*cos(theta(3)) + sin(theta(1))*cos(theta(2))*sin(theta(3))

  CBI(1,2) = -cos(theta(2))*sin(theta(3)) + sin(theta(1))*sin(theta(2))*cos(theta(3))
  CBI(2,2) =  cos(theta(1))*cos(theta(3))
  CBI(3,2) =  sin(theta(2))*sin(theta(3)) + sin(theta(1))*cos(theta(2))*cos(theta(3))

  CBI(1,3) =  cos(theta(1))*sin(theta(2))
  CBI(2,3) = -sin(theta(1))
  CBI(3,3) =  cos(theta(1))*cos(theta(2))

end function CBI

! transformation matrix between angular velocity between intertial and body axis
function SIB(theta)
  real(dp), intent(in)    :: theta(num_spat_dim)
  real(dp)                :: SIB(num_spat_dim, num_spat_dim)

  SIB(1,1) = cos(theta(3))
  SIB(2,1) = cos(theta(1))*sin(theta(3))
  SIB(3,1) = 0.0_dp

  SIB(2,1) = -sin(theta(3))
  SIB(2,2) = cos(theta(1))*cos(theta(3))
  SIB(3,2) = 0.0_dp

  SIB(3,1) = 0.0_dp
  SIB(3,2) = -sin(theta(1)) 
  SIB(3,3) = 1.0_dp

end function SIB

! function to create body and set the class variables
!!$function makeBody()
!!$type(body) :: makebody
!!$end function makeBody
  
end module body_mod
