!===========================================================
! program to run dynamics simulation of a physical pendulum
!===========================================================
program pendulum

  use constants
  use utils
  use body_mod
  use matrix_utils

  implicit none

!  type(body) :: alpha
  real(dp)   :: residual(12)
  real(dp)   :: jacob(12,12)

  ! solver related variables
  integer neq,info(15),idid,lrw,liw,iwork(1000),ipar, my, ii
  double precision t, y(12), yprime(12), tout, rtol(12),atol(12), rwork(1055), rpar, h0
  external RES , JAC

  write(*,*) "Setting up the test problem"

  call create_body()

  print*, alpha

  stop

  ! set number of equations to solve
  NEQ=12

  IDID=0  !?

  ! set tolerances
  Rtol = 0.005
  Atol = 0.005

  ! intial time
  t=0.00

  ! final time
  tout=0.03


  ! initial state

  ! info block
  do ii=1,15
     info(ii)=0    
  end do
  info(2)=1
  info(3)=1

  ! allocate size of work arrays
  liw=1000
  lrw=1055 


  ! (1) create a pendulum body

  ! (2) set the res and jacobian for the linear system

  ! (3) call the solver
  write(*,*) '***************'
  write(*,*) 'time,                              positions     ' 
  do  while (t .le. tout) 
     call ddassl(res,neq,t,y,yprime,tout,info,rtol,atol, idid,rwork,lrw,iwork,liw,rpar,ipar,jac)
     write(*,'(8f13.8)') t,y(1),y(2),y(3),y(4),y(5),y(6),y(7)
  end do

end program pendulum

!********************************
! implementation of the residual
!********************************
subroutine res(t,y,yprime,delta,ires,rpar,ipar)
  ! implements the residual of the governing equation at the for the given t, y, yprime values
  use constants
  use utils,only:get_vector
  use body_mod,only:alpha, R_rigid
!  use jacobian_mod, only: res_rigid
  implicit none

  real(dp)  :: t
  real(dp)  :: y(12)
  real(dp)  :: yprime(12)
  real(dp), intent(out)  :: delta(12)
  integer(sp) :: ires
  real(dp)  :: rpar(*)
  integer(sp) :: ipar(*)

  integer(sp):: i,j,k

  
  ! sanity check 
!!$  do i = 1, 12
!!$     if (y(i) .lt. dzero)  ires = -1; return
!!$  end do

  ! option to compute yprime by itself

  !?? get the residual
  delta = get_vector(R_rigid(alpha),4)
  print*, "residual vector = ", delta

end subroutine res

!********************************
! implementation of the jacobian
!********************************
subroutine jac
  ! JACOBIAN  = jac_rigid(alpha,a) !?? implement jacobian
  stop"where is the impl"
end subroutine jac






!  frame_A%rot_mat     = rotmat(frame_A%theta_alpha)
!  body_A%v_alpha   = matmul(frame_A%rot_mat, frame_A%v_alpha)

!  print*, frame_A%v_alpha
!  print*, body_A%v_alpha
!  print*, ""

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

subroutine create_body()

  use constants
  use utils
  use body_mod
  use matrix_utils

  implicit none

  real(dp), parameter           :: m  = 1.0_dp
  type(vector), parameter       :: g0 = vector((/ dzero, dzero, -1.0_dp /))
  type(vector), parameter       :: re = vector((/ 1.0_dp, 1.0_dp, 1.0_dp /)) ! position vector of a point on the circumference

  real(dp)                      :: r(num_spat_dim), theta(num_spat_dim), v(num_spat_dim), omega(num_spat_dim)

!!$
!!$!(r, theta, v, omega, r_dot, theta_dot, v_dot, omega_dot,&
!!$!     & mass, c, J, p, h, K, M, C_mat, S, S_dot, fr, gr)
!!$

  ! define positions of body frame from inertial (q terms)
  r       = (/ 2.0_dp, 2._dp, 2.0_dp /)
  theta   = (/ deg2rad(30.0d0), deg2rad(45.0d0), deg2rad(60.0d0) /)
  v     = (/ 0.0d0, 0.0d0, 0.0d0 /)
  omega = (/ 0.0d0, 0.0d0, 1.0d0 /)

  alpha%r         = vector(r)
  alpha%theta     = vector(theta)
  alpha%v         = vector(v)
  alpha%omega     = vector(omega)

  ! define velocities of the body frame (q_dot terms)
  alpha%r_dot     = zeroV
  alpha%theta_dot = zeroV
  alpha%v_dot     = zeroV
  alpha%omega_dot = zeroV

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

end subroutine create_body
