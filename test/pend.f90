program pendulum
! program to run dynamics simulation of a pendulu use constants
  use utils
  use body_mod
  implicit none

  type(body) :: alpha
  real(dp)   :: residual(12)
  real(dp)   :: jac(12,12)
  
  ! (1) create a pendulum body
  
  
  ! (2) Set the res and jacobian for the linear system



  ! (3) Call the solver

  
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
!!$
!!$SUBROUTINE RES(T,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR, senpar)
!!$
!!$END SUBROUTINE RES
