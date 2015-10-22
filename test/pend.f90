!===========================================================
! program to run dynamics simulation of a physical pendulum
!===========================================================
program pendulum
  use dispmodule 
  use constants
  use utils
  use body_mod
  use matrix_utils
  use jacobian_mod

  implicit none
  real(dp)         :: delta(12)
  real(dp)         :: pd(12,12)
  ! solver related variables
  integer(sp)      :: neq,info(15),idid,lrw,liw,iwork(1000),ipar, my, ii
  real(dp)         :: t, y(12), yprime(12), tout, rtol(12),atol(12), rwork(1055), rpar
  ! state variables
  real(dp)         :: r(num_spat_dim), theta(num_spat_dim), v(num_spat_dim), omega(num_spat_dim) ! state vars
  real(dp)         :: r_dot(num_spat_dim), theta_dot(num_spat_dim), v_dot(num_spat_dim), omega_dot(num_spat_dim) ! state vars
  ! other variables
  real(dp)         :: m, g0(num_spat_dim), re(num_spat_dim)

  external RES , JAC

  write(*,*) "Setting up the test problem"

  ! ******************************************************
  ! (1) create a pendulum body
  ! ******************************************************

  r       = (/ 1.0_dp, 2._dp, 3.0_dp /)
  theta   = (/ deg2rad(10.0d0), deg2rad(-20.0d0), deg2rad(-40.0d0) /)
  v       = (/ 1.0d0, -2.0d0, 3.0d0 /)
  omega   = (/ -1.0d0, 3.0d0, -4.0d0 /)

  r_dot       = (/ 1.0d0, -2.0d0, 3.0d0 /)
  theta_dot   = (/ -1.0d0, 2.0d0, -1.0d0 /)
  v_dot       = (/ 1.0d0, 2.0d0, -2.0d0 /)
  omega_dot   = (/ -1.0d0, 4.0d0, -2.0d0 /)

  ! initial state values
  Y       = (/r, theta, v, omega /)
  YPRIME  = (/r_dot, theta_dot, v_dot, omega_dot /)

  !  call disp('   Y      =  ', Y)
  !  call disp('   Yprime =  ',yprime)

  m  = 5.0d0

  g0 = (/ -1.0d0, -1.0d0, -1.0d0/)
  re = (/ -1.0d0, 2.0d0, 1.0d0/)

  !  CALL DISP('   g0 =   ', g0, SEP=', ', ORIENT = 'ROW') 
  !  CALL DISP('   re =   ', re, SEP=', ', ORIENT = 'ROW') 
  
  call create_body(m, vector(g0), vector(re), Y, YPRIME, alpha)

!!$  delta = get_vector_elements(R_rigid(alpha),4)
!!$  call print_body(alpha)
!!$  call disp("   R   =   ", delta)
!!$
!!$  pd = jac_rigid (alpha, 2.0d0)
!!$  call disp("   JAC =   ", pd)
!!$
!!$  stop

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
  tout=1.0d-1

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

     !     write(*,'(8f13.2)') t, y(1), y(2), y(3), y(4), y(5), y(6), y(7), y(8), y(9), y(10), y(11), y(12)

     call disp('      t, y(1:12) = ',(/ t, y/), SEP=', ', ORIENT = 'ROW')

  end do

end program pendulum

!********************************
! implementation of the residual
!********************************
subroutine res(t,y,yprime,delta,ires,rpar,ipar)
  ! implements the residual of the governing equation at the for the given t, y, yprime values
  use dispmodule 
  use constants
  use utils,only:get_vector_elements,disp
  use body_mod
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
  call SetStateVars(Y, YPRIME, alpha)
  !?? get the residual
  delta = get_vector_elements(R_rigid(alpha),4)
  !  call print_body(alpha)
  call disp("   R   =   ", delta)
  !  stop
  !  call print_body(alpha)
  !  print*,"body=",alpha
  !  print*, "residual vector = ", delta
  !  call disp(t,delta)

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
  use dispmodule 
  use constants,only:dp,sp
  use utils,only:disp_mat
  use body_mod
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

  !  print*,"a=", cj
  !  cj = 1.0d-1
  call SetStateVars(Y, YPRIME, alpha)
  pd = jac_rigid (alpha, cj)
  call disp("   JAC =   ", pd)
  !  print*, "body =" , alpha
  !  call print_body(alpha)
  ! print*, "jacobian matrix ="
  ! call disp_mat(t, pd, 12, 12)
  ! stop"where is the impl"
  return
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
