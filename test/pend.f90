program pendulum
  ! program to run dynamics simulation of a pendulu use constants
  use constants
  use utils
  use body_mod
  use matrix_utils

  implicit none

  type(body) :: alpha
  real(dp)   :: residual(12)
  real(dp)   :: jacob(12,12)
  INTEGER NEQ,INFO(15),IDID,LRW,LIW,IWORK(1000),IPAR, MY, II
  DOUBLE PRECISION T, Y(12), YPRIME(12), TOUT, RTOL(12),ATOL(12), RWORK(1055), RPAR,H0
  double precision :: temp

  real(dp), parameter           :: m = done
  type(vector), parameter       :: g0 = vector((/ dzero, dzero, -1.0_dp /))

  type(vector), parameter       :: zeroV = vector((/ dzero, dzero, dzero /))
  type(vector), parameter       :: unitV = vector((/ dzero, dzero, dzero /))

  type(matrix)                  :: idM   
  type(matrix)                  :: zeroM 
  type(matrix)                  :: unitM 

  type(vector), parameter       :: re = vector((/ 1.0_dp, 1.0_dp, 1.0_dp /)) ! position vector of a point on the circumference of the sphere
!  type(frame)                   :: frame_A

  external RES , JAC
 
!  frame_A%rot_mat     = rotmat(frame_A%theta_alpha)
!  body_A%v_alpha   = matmul(frame_A%rot_mat, frame_A%v_alpha)

!  print*, frame_A%v_alpha
!!  print*, body_A%v_alpha
!  print*, ""

  idM   = matrix(eye(num_spat_dim))
  zeroM = matrix(zeros(num_spat_dim))
  unitM = matrix(ones(num_spat_dim))


  ! define positions of body frame from inertial
  alpha%r         = vector((/ 2.0_dp, 2._dp, 2.0d0 /))
  alpha%theta     = vector((/ deg2rad(30.0d0), deg2rad(45.0d0), deg2rad(60.0d0)/))

  ! define velocities of the body frame
  alpha%r_dot     = zeroV
  alpha%omega_dot = zeroV

  ! rotation and rate matrices


  ! mass of the body
  alpha%m  = m

  ! mass moment of intertia (second moment)
  temp     = (2.0_dp*m*abs(re)**2)/5.0_dp  
  alpha%J  = temp*idM
  
  ! first moment of mass
  alpha%c  = m*re
  
  ! force 
  alpha%fr = m*g0
  
  ! torque
  alpha%gr = alpha%J*alpha%omega_dot

  ! set number of equations to solve
  NEQ=12

  IDID=0  !?

  ! set tolerances
  Rtol = 0.005
  Atol = 0.005

  ! intial time
  T=0.00

  ! final time
  TOUT=0.03


  ! initial state

  ! info block
  DO II=1,15
     INFO(II)=0    
  end do
  INFO(2)=1
  INFO(3)=1

  ! allocate size of work arrays
  LIW=1000
  LRW=1055 


  ! (1) create a pendulum body

  ! (2) Set the res and jacobian for the linear system

  ! (3) Call the solver
  WRITE(*,*) '***************'
  WRITE(*,*) 'Time,                              Positions     ' 
  do  while (t .le. tout) 
     CALL DDASSL(RES,NEQ,T,Y,YPRIME,TOUT,INFO,RTOL,ATOL, IDID,RWORK,LRW,IWORK,LIW,RPAR,IPAR,JAC)
     WRITE(*,'(8F13.8)') T,Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7)
  end do


!!$
!!$  real(dp), parameter       :: t0=0, t_final=1
!!$  real(dp)                  :: current_time
!!$  integer, parameter        :: nSteps=10 ! number of time steps
!!$  integer                   :: i,j,k
!!$  real(dp)                  :: q0(dim_q)
!!$  real(dp)                  :: q(0:nSteps,dim_q)
!!$  real(dp)                  :: q_dot(0:nSteps,dim_q)
!!$  real(dp)                  :: q_double_dot(0:nSteps, dim_q)  ! state variable, its first and second derivatives
!!$  real(dp), parameter       :: tol = 1.0d-8
!!$  real(dp)                  :: residual
!!$  real(dp)                  :: del_q(dim_q)
!!$  
!!$  type(body)                :: body_A
!!$  type(body_fixed_frame)    :: frame_A
!!$  type(matrix) :: test(1,1)
!!$  type(vector) ::dd(2)



end program pendulum


SUBROUTINE RES(T,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
  ! implements the residual of the governing equation at the for the given T, Y, YPRIME values
  use constants
  use utils
  implicit none

  real(dp)  :: t
  real(dp)  :: y(12)
  real(dp)  :: yprime(12)
  real(dp), intent(out)  :: delta(12)
  integer(sp) :: IRES
  real(dp)  :: RPAR(*)
  integer(sp) :: IPAR(*)

  integer(sp):: i,j,k

!!$  do i = 1, 12
!!$     if (y(i) .lt. dzero)  ires = -1; return
!!$  end do

  ! option to compute yprime by itself


END SUBROUTINE RES

! implement th jacoain if needed
subroutine jac
  stop"where is the impl"
end subroutine jac

! rotates a vector from body frame inertial
function CBI(theta)
  use constants
  real(dp), intent(in)    :: theta(num_spat_dim)
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
  use constants
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
