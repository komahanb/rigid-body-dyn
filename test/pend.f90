!===========================================================
! program to run dynamics simulation of a physical pendulum
!===========================================================
program pendulum
  ! import the necessary modules
  !  use iso_c_binding
  use global_constants
  use global_variables
  use types
  use utils
  use bodies, only: create_body
!  use joints, only: joint
  use residual, only: get_residual
  use jacobian, only: get_jacobian
  
  implicit none

  real(dp), dimension(NUM_SPAT_DIM)      :: re    ! position of the point of interest in the body in the body frame
  real(dp)                               :: mass ! mass
  real(dp)                               :: qr(NUM_STATES)     = 0.0_dp
  real(dp)                               :: qr_dot(NUM_STATES) = 0.0_dp
  real(dp)                               :: res(NUM_STATES), jac(NUM_STATES, NUM_STATES)
  real(dp)                               :: dq = 0.0_dp
  type(body)                             :: body1

  call disp("==================================")
  call disp("-------Rigid body dynamics--------")
  call disp("==================================")

  ! ******************************************************
  ! (1) set the initial states and attributes of the body
  ! ******************************************************

  call disp(" >> Setting up the test problem...")

  ! define the initial states
  qr(1:3)    = (/ 1.0d0, 1.0d0, 1.0d0 /)
  qr(4:6)    = (/ deg2rad(20.0d0), deg2rad(22.0d0), deg2rad(45.0d0) /)
  qr(7:9)    = (/ 1.0d0,  2.0d0, 1.0d0 /)
  qr(10:12)  = (/ 2.0d0, 4.0d0, 5.0d0 /)

  ! define the intial time derivative of state
  qr_dot(1:3)   = (/ 2.0d0, 2.0d0, 2.0d0 /)
  qr_dot(4:6)   = (/ 5.0d0, 4.0d0, 2.0d0 /)
  qr_dot(7:9)   = (/ 0.0d0, 0.0d0, 0.0d0 /)
  qr_dot(10:12) = (/ 0.0d0, 0.0d0, 0.0d0 /)

  ! define the attributes of the body
  mass      = 2.0d0
  re        = (/ 1.0d0, 1.0d0, 1.0d0 /)  

  ! ******************************************************
  ! (2) create a pendulum body using the state and attrs
  ! ******************************************************

  call disp(" >> Creating a body...")
  body1 = create_body(mass, vector(re), qr, qr_dot)
  call disp(" >> Body created successfully...")

  !  call testtotalderivative(5, 1.0d-6)
  
!!$  ! ******************************************************
!!$  ! (3) Assemble the residual and jacobian using the body
!!$  ! ******************************************************
!!$
!!$  call disp(" >> Assembling the residuals...")
!!$  res  = get_residual(body1)
!!$  call disp(" >> Residual assembly complete...")
!!$  
!!$  call disp("   R   =   ", res)
!!$
!!$  call disp(" >> Assembling the Jacobian...")
!!$  jac  = get_jacobian(body1)
!!$  call disp(" >> Jacobian assembly complete...")
!!$
!!$  call disp("   JAC =   ",jac)
!!$
!!$  ! ******************************************************
!!$  ! (3) Solve the linear system and compute update delta_q
!!$  ! ******************************************************
!!$
!!$  call disp(" >> Calling the linear solver...") 
!!$!  dq = linear_solve(jac,-res,'GMRES')
!!$  call disp(" >> Update to residual is computed...") 
!!$
!!$  ! ******************************************************
!!$  ! (4) Update the state variables and then update the body
!!$  ! ******************************************************

  


end program pendulum
