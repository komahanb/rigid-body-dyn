!=====================================================================
! program to run test the self consistency of the linear system
!=====================================================================
program test_linsys

  ! import the necessary modules
  !  use iso_c_binding ! to use tripan
  use global_constants
  use global_variables
  use types
  use utils
  use solver_utils
  use bodies, only: create_body, set_state
  use residual, only: get_residual
  use jacobian, only: get_jacobian
  use finite_diff
  implicit none

  ! position of the point of interest on the body measured in body frame
  real(dp), dimension(NUM_SPAT_DIM)  :: re    
  ! mass of the body
  real(dp)                           :: mass 
  ! the state
  real(dp)                           :: q(TOT_NDOF)= 0.0_dp
  real(dp)                           :: qtmp(TOT_NDOF)= 0.0_dp

  ! time derivative of state
  real(dp)                           :: q_dot(TOT_NDOF) = 0.0_dp
  real(dp)                           :: q_dot_tmp(TOT_NDOF) = 0.0_dp

  ! residual
  real(dp)                           :: res(TOT_NDOF)

  real(dp)                           :: res2(TOT_NDOF)

  ! jacobian
  real(dp)                           :: jac(TOT_NDOF, TOT_NDOF)
  real(dp)                           :: jac2(TOT_NDOF, TOT_NDOF)

  ! update from newton iterations
  real(dp)                           :: dq(TOT_NDOF) = 0.0_dp
  ! body -- single pendulum
  type(body)                         :: body1

  ! loop variable
  integer(sp)                        :: k

  call disp("========================================================")
  call disp("                 'TEST jacobain '                       ")
  call disp("========================================================")

  ! state vector
  CALL RANDOM_NUMBER(q);  CALL RANDOM_NUMBER(q_dot) ;  
  q_dot = q_dot*q_dot + 1.0_dp;
  
  ! inertial props
  mass  = 2.0_dp;   re   = (/ 0.1_dp, 0.2_dp, 0.3_dp /);
  body1 = create_body(mass, vector(re), q, q_dot);

  ! co-eff for jacobian
  aa = 1.0d0 
  
  jac = get_jacobian(body1) !actual
  jac2 = finite_difference2(q, q_dot, aa, 1.0d-6) !finite diff
  

  call disp("exact=", jac) 
  print*,""
  call disp("appro=", jac2)
  print*,""
  call disp("DIFF=", jac2-jac)
  print*,""
  call disp("Max Diff=", maxval(jac2-jac))

end program test_linsys

! impl for finite-differencing 
subroutine residual2(q, qdot, f)

  use global_constants, only: sp, dp, TOT_NDOF, NUM_SPAT_DIM
  use types
  use bodies
  use residual

  real(dp), intent(in), dimension(TOT_NDOF)  :: q, qdot
  real(dp), intent(out), dimension(TOT_NDOF) :: f
  real(dp), dimension(NUM_SPAT_DIM)  :: re    
  real(dp)   :: mass
  type(body) :: body1

  mass  = 2.0_dp;   re   = (/ 0.1_dp, 0.2_dp, 0.3_dp /);
  body1 = create_body(mass, vector(re), q, qdot);

  f = get_residual(body1);  

end subroutine residual2

subroutine residual()
end subroutine residual
