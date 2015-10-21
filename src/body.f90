! module that defines the configuration of a single body
module body_mod
  use utils
  use matrix_utils
  implicit none

  private
  
  public body, CBI, SIB
  
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

contains

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
  
end module body_mod
