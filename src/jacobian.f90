! moduele that assembles the jacobian
! terms in jacobian come from (a) body (b) joints (c)
module jacobian_mod
use utils
use body_mod
implicit none

type jacobian
   type(matrix) :: D_R (4,4)
   type(matrix) :: S_R (4,4)
   type(matrix) :: S_S(1,1)
   type(matrix) :: S_RS(4,1)
   type(matrix) :: S_SR(1,4)
end type jacobian

contains

  ! decompose the jacobian in to double precision matrix
 function getJacobian(jac)
   type(jacobian)  :: Jac
   real(dp)        :: getJacobian(12,12)
   type(matrix)   :: D_R (4,4)
   D_R = Jac%D_R
   
   getJacobian = getJacobian_2d(D_R,4,4)

   print*, decompose

 end function getJacobian
!!$ 
!!$ function residual(res)
!!$   real(dp) :: residual(size(res),num_spat_dim)
!!$ end function residual
!!$

!!$function Jac(alpha)
!!$  type(body)   :: alpha
!!$  type(matrix) :: Jac(2,2) ! where 4 is the state vector length
!!$
!!$!  Jac(1,1) =  D_R(alpha) + S_R(alpha) ! 4x4 ! 12*12
!!$  Jac(2,1) =  S_RS(alpha)
!!$
!!$!  Jac(1,2) =  S_SR(alpha)
!!$!  Jac(2,2) =  S_S(alpha)
!!$  
!!$end function Jac


! find the jacobian of the equation of motion for the supplied BODY alpha
function D_R(alpha, a)

  type(body)    :: alpha
  real(dp)      :: a
  type(matrix)  :: D_R(4,4) ! where 4 is the state vector length
  type(matrix)  :: O ! zero matrix
  type(matrix)  :: U ! unit matrix
  type(matrix)  :: I ! identity matrix

  ! some useful matrices
  O = matrix(zeros(num_spat_dim))
  U = matrix(ones(num_spat_dim))
  I = matrix(eye(num_spat_dim))
  
  ! assemble jacobian
  D_R(1,1) = a * alpha%C_mat
  D_R(2,1) = skew(alpha%C_mat * alpha%r_dot) * alpha%S
  D_R(3,1) = -1.0_dp*U
  D_R(4,1) = O

  D_R(1,2) = O
  D_R(2,2) = alpha%S_dot + skew(alpha%S*alpha%theta_dot)*alpha%S + a * alpha%S
  D_R(3,2) = O
  D_R(4,2) = -1.0_dp*U

  D_R(1,3) = O
  D_R(2,3) = O
  D_R(3,3) = alpha%m*(a*U + skew(alpha%omega))
  D_R(4,3) = -a*skew(alpha%c) + skew(skew(alpha%c)*alpha%omega)  - alpha%m*skew(alpha%V) - skew(alpha%omega)*skew(alpha%C)

  D_R(1,4) = O
  D_R(2,4) = O
  D_R(3,4) = a*skew(alpha%C) + skew(alpha%C) * skew(alpha%omega)
  D_R(4,4) = a*alpha%J  + skew(alpha%omega)*alpha%J - skew(alpha%J*alpha%omega) -skew(alpha%c)*skew(alpha%V)

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
function S_R(alpha, a)

  type(body)   :: alpha
  real(dp)     :: a
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

! linear combination of stiffness and mass matrix
function S_S(alpha, a)
  type(body)   :: alpha
  real(dp)     :: a
  type(matrix) :: S_S

  S_S = alpha%K  + alpha%b*alpha%M ! K + b M

end function S_S

function S_RS(alpha, a)
  type(body)   :: alpha
  real(dp)     :: a
  type(matrix) :: S_RS(4,1)
  type(matrix) :: O ! zero matrix

  O = matrix(zeros(num_spat_dim))

  S_RS(1,1) = O
  S_RS(2,1) = O
  S_RS(3,1) = alpha%b*alpha%p + a*skew(alpha%omega)*alpha%p !bp +a wx p
  S_RS(4,1) = alpha%b*alpha%h + a*(skew(alpha%v)*alpha%p  + skew(alpha%omega)*alpha%h) !bh +a(vx p + wx h)

end function S_RS


function S_SR(alpha, a)
  type(body)   :: alpha
  real(dp)     :: a
  type(matrix) :: S_SR(1,4)
  type(matrix) :: O ! zero matrix

  O = matrix(zeros(num_spat_dim))

  S_SR(1,1) = O
  S_SR(1,2) = O
  S_SR(1,3) = a * trans(alpha%p) !a p^T
  S_SR(1,4) = a * trans(alpha%p) !a h^T
 
end function S_SR


end module jacobian_mod
