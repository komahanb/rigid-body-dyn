program pendulum
  ! program to run dynamics simulation of a pendulu use constants
  use utils
  use body_mod
  implicit none

  type(body) :: alpha
  real(dp)   :: residual(12)
  real(dp)   :: jacob(12,12)

  INTEGER NEQ,INFO(15),IDID,LRW,LIW,IWORK(1000),IPAR, MY, I, II
  DOUBLE PRECISION T, Y(20), YPRIME(20), TOUT, RTOL(12),ATOL(12), RWORK(1055), RPAR,H0

  external RES , JAC

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
  DO I=1,15
     INFO(I)=0    
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
