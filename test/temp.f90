!===========================================================
! program to run dynamics simulation of a physical pendulum
!===========================================================
program pendulum

  ! import the necessary modules
  use global_constants
  use global_variables
  use types
  
  implicit none

  real(dp), dimension(NUM_SPAT_DIM)      :: re    ! position of the point of interest in the body in the body frame
  real(dp)                               :: mass ! mass

  real(dp), dimension(NUM_STATES)        :: qr(NUM_STATES) = 0.0_dp    
  real(dp), dimension(NUM_STATES)        :: qr_dot(NUM_STATES) = 0.0_dp    

  !type(body), dimension(MAX_NUM_BODIES)  :: body
  type(body)                              :: body

  call disp("==================================")
  call disp("-------Rigid body dynamics--------")
  call disp("==================================")

  ! ******************************************************
  ! (1) set the initial states and attributes of the body
  ! ******************************************************

  call disp(" >> Setting up the test problem...")

  ! define the initial states
  qr(1:3)    = (/ 0.0d0, 0.0d0, 0.0d0 /)
  qr(4:6)    = (/ deg2rad(0.0d0), deg2rad(0.0d0), deg2rad(45.0d0) /)
  qr(7:9)    = (/ 0.0d0,  0.0d0, 0.0d0 /)
  qr(10:12)  = (/ 0.0d0, 0.0d0, 0.0d0 /)

  ! define the intial time derivative of state
  qr_dot(1:3)   = (/ 0.0d0, 0.0d0, 0.0d0 /)
  qr_dot(4:6)   = (/ 0.0d0, 0.0d0, 0.0d0 /)
  qr_dot(7:9)   = (/ 0.0d0, 0.0d0, 0.0d0 /)
  qr_dot(10:12) = (/ 0.0d0, 0.0d0, 0.0d0 /)

  ! define the attributes of the body
  mass      = 1.0d0
  re        = (/ 0.0d0, 0.0d0, 0.0d0  /)  

  ! ******************************************************
  ! (2) create a pendulum body using the state and attrs
  ! ******************************************************

  call disp(" >> Creating a body...")

  !  call create_body(m, vector(g), vector(r), q, qr, alpha)

  ! we can just create an array of bodies for multiple bodies
  body1 = create_body(mass, vector(re), qr, qr_dot)

  ! ******************************************************
  ! (3) Assemble the residual and jacobian using the body
  ! ******************************************************

  call disp(" >> Assembling the residuals...")

  res  = assembleResidual(body1)

  call disp(" >> Residual assembly complete...")
  !  call disp("   R   =   ", res)

  call disp(" >> Assembling the Jacobian...")

  jac  = assembleJacobian(body1)

  call disp(" >> Jacobian assembly complete...")
  !  call disp("   JAC =   ",jac)

  ! ******************************************************
  ! (3) Solve the linear system and compute update delta_q
  ! ******************************************************

  call disp(" >> Solving the linear system...")

  dq = linear_solve(jac,-res,'GMRES')

  ! ******************************************************
  ! (4) Update the state variables and then update the body
  ! ******************************************************

  ! update the state
  q = q + dq;
  qdot = qdot + a *dq;   

  call update_body_state(q, qr, body1)


  print*,"dq=", dq

  stop

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
  info(5) = 0   ! Yes - Set INFO(5)=0   ! No  - Set INFO(5)=1 ! provide JAC

  info(11) = 0 ! 1 = Y and Y prime are consistent ; 0 = not consistent compute automatically

  ! allocate size of work arrays
  liw=1000;  lrw=1055 ;

  ! (3) call the solver
  do while (t .le. tout) 


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
  call update_state_vars(Y, YPRIME, alpha)

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

  call SetStateVars(Y, YPRIME, alpha)
  pd = jac_rigid (alpha, cj)
  call disp("   JAC =   ", pd)

  return
end subroutine jac
