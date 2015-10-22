!===========================================================
! program to run dynamics simulation of a physical pendulum
!===========================================================
program pendulum

  use constants
  use utils
  use body_mod
  use matrix_utils

  implicit none

  ! solver related variables
  integer neq,info(15),idid,lrw,liw,iwork(1000),ipar, my, ii
  double precision t, y(12), yprime(12), tout, rtol(12),atol(12), rwork(1055), rpar, h0

  external RES , JAC

  write(*,*) "Setting up the test problem"

  ! ******************************************************
  ! (1) create a pendulum body
  ! ******************************************************

  call create_body()
  
  ! ******************************************************
  ! (2) set the res and jacobian for the linear system
  ! ******************************************************

  ! set number of equations to solve
  NEQ=12

  IDID=0  

  ! set tolerances
  Rtol = 0.05
  Atol = 0.05

  ! intial time
  t=0.00

  ! final time
  tout=0.03

  ! info block
  do ii=1,15
     info(ii)=0    
  end do
  info(2)=1
  info(3)=1

  ! derivatives automatically by numerical differences
  info(5) = 1   ! Yes - Set INFO(5)=0   ! No  - Set INFO(5)=1 ! provide JAC

  info(11) = 0 ! 1 = Y and Y prime are consistent ; 0 = not consistent compute automatically

  ! allocate size of work arrays
  liw=1000;  lrw=1055 ;

  ! (3) call the solver
  write(*,*) '************************************************************************************&
       &*************************************************************************************'
  write(*,*) '       Time', '            q1', '           q2', '          q3', '            q4', &
       &'          q5', '            q6', '           q7', '          q8',&
       & '          q9', '           q10', '         q11',  '           q12'
  write(*,*) '************************************************************************************&
       &*************************************************************************************'
  do while (t .le. tout) 
     call ddassl(res,neq,t,y,yprime,tout,info,rtol,atol,&
          & idid,rwork,lrw,iwork,liw,rpar,ipar,jac)
     write(*,'(8f13.2)') t,y(1),y(2),y(3),y(4),y(5),y(6),y(7),y(8),y(9),y(10),y(11),y(12)
  end do

end program pendulum

!********************************
! implementation of the residual
!********************************
subroutine res(t,y,yprime,delta,ires,rpar,ipar)
  ! implements the residual of the governing equation at the for the given t, y, yprime values
  use constants
  use utils,only:get_vector_elements
  use body_mod,only:alpha, R_rigid
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
  delta = get_vector_elements(R_rigid(alpha),4)
!  print*, "residual vector = ", delta

  call disp(t,delta)

end subroutine res

!********************************
! implementation of the jacobian
!********************************
subroutine jac(t,y,yprime,pd,cj,rpar,ipar)
!!$C  to define the matrix of partial derivatives
!!$C             PD=DG/DY+CJ*DG/DYPRIME
!!$C         CJ is a scalar which is input to JAC.
!!$C         For the given values of T,Y,YPRIME, the
!!$C         subroutine must evaluate the non-zero partial
!!$C         derivatives for each equation and each solution
!!$C         component, and store these values in the
!!$C         matrix PD.  The elements of PD are set to zero
!!$C         before each call to JAC so only non-zero elements
!!$C         need to be defined.
!!$C
!!$C         Subroutine JAC must not alter T,Y,(*),YPRIME(*), or CJ.
!!$C         You must declare the name JAC in an EXTERNAL statement in
!!$C         your program that calls DDASSL.  You must dimension Y,
!!$C         YPRIME and PD in JAC.

  use constants,only:dp,sp
  use body_mod,only:alpha
  use jacobian_mod,only:jac_rigid
  implicit none
  real(dp)     :: cj
  real(dp)     :: t
  real(dp)     :: y(12)
  real(dp)     :: yprime(12)
  integer(sp)  :: ires
  real(dp)     :: rpar(*)
  integer(sp)  :: ipar(*)
  real(dp)     :: pd(12,12)
  integer(sp)  :: i,j,k

  ! JACOBIAN  = jac_rigid(alpha,a) !?? implement jacobian
 ! print*,cj 
!  stop
  cj = 0.1
  pd = jac_rigid (alpha, cj)

  print*,  shape(pd)
  
  call disp_mat(t, pd, 12, 12)

  stop

  !  stop"where is the impl"

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
  r       = (/ 0.0_dp, 0._dp, 0.0_dp /)
  theta   = (/ deg2rad(10.0d0), deg2rad(0.0d0), deg2rad(0.0d0) /)
  v     = (/ 1.0d0, 0.0d0, 0.0d0 /)
  omega = (/ 1.0d0, 0.0d0, 0.0d0 /)

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



subroutine disp(t,y)

  use constants,only:dp,sp
  implicit none

  real(dp)    :: t, y(12)
  integer(sp) :: i

  write(*,'(13f13.2)') t, (y(i),i=1,12)


end subroutine disp


subroutine disp_mat(t, A,m,n)

  use constants, only:dp,sp
  implicit none

  real(dp), intent(in) :: t
  integer(sp)          :: i, j, m, n
  real(dp)             :: A(m,n)


  do i = 1, n
     write(*,'(13f13.2)') t, (A(j,i),j=1,m)
  end do


end subroutine disp_mat
