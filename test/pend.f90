!===========================================================
! program to run dynamics simulation of a physical pendulum
!===========================================================
program pendulum

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
  ! time, and norms to stop integration
  real(dp)                           :: time,l2_norm, linf_norm, dh
  ! loop variable
  integer(sp)                        :: k


  call disp("========================================================")
  call disp("                 Rigid body dynamics                    ")
  call disp("========================================================")

  ! init method that has sanity check on the input settings
  ! has nothing as of now, we can use to load input files, python etc
  call init() 
  
  !-------------------------------------------------------------------!
  ! (1) set the initial states and attributes of the body
  !-------------------------------------------------------------------!

  call disp(" >> Setting up the test problem...")
  
!!$  ! define the initial states
!!$  q(1:3)    = (/ 1.0_dp, 2.0_dp, 3.0_dp /)
!!$  q(4:6)    = (/ deg2rad(20.0_dp), deg2rad(0.0_dp), deg2rad(20.0_dp) /)
!!$  q(7:9)    = (/ 1.0_dp,  2.0_dp, 1.0_dp /)
!!$  q(10:12)  = (/ 2.0_dp, 4.0_dp, 5.0_dp /)
!!$
!!$  ! define the intial time derivative of state
!!$  q_dot(1:3)   = (/ 2.0_dp, 2.0_dp, 2.0_dp /)
!!$  q_dot(4:6)   = (/ 0.25_dp, 0.0_dp, 0.0_dp /)
!!$  q_dot(7:9)   = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
!!$  q_dot(10:12) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)

!  CALL init_random_seed()         ! see example of RANDOM_SEED
  CALL RANDOM_NUMBER(q)
  CALL RANDOM_NUMBER(q_dot)  
  q_dot = q_dot*q_dot+1.0_dp

  ! define the attributes of the body
  mass      = 2.0_dp
  ! point of interest on the body  is the point where the body axis is located
  CALL RANDOM_NUMBER(re)
  !  re        = (/ 0.0_dp, 0.0_dp, 0.0_dp /)   

  !-------------------------------------------------------------------!
  ! (2) create a pendulum body using the state and attr above
  !-------------------------------------------------------------------!
  aa = 1.0_dp

  call disp(" >> Creating a body...")

  body1 = create_body(mass, vector(re), q, q_dot)

  call disp(" >> Body created successfully...")
    
  ! residual assembly
  call disp(" >> Assembling the residuals...")
  res  = get_residual(body1)
  call disp(" >> Residual assembly complete...")
  call disp("   R   =   ", res)

  ! jacobian assembly
  call disp(" >> Assembling the Jacobian...")
  jac  = get_jacobian(body1)
  call disp(" >> Jacobian assembly complete...")
  call disp("   JAC =   ",jac)


  ! implement FD to check if the jacobian is the same
  qtmp = q
  call set_state(q, q_dot, body1)
  dh = 1.0d-8

  res = get_residual(body1)
     
  var: do  k = 1 , TOT_NDOF

     ! perturb k-th variable
     qtmp( k ) = qtmp (k)  + dh
     ! updating the body state
     call set_state(qtmp, q_dot, body1)
     ! find the solution
     res2 = get_residual(body1)
     !    call disp("   R   =   ", res)
     
     dq =  (res2-res)/dh
     
     jac2 (:, k) = dq
     
     qtmp = q

  end do var
  call disp("   dJAC =   ",jac2)
 
  call disp("   Error =   ",jac2-Jac)



  stop"stopped"

  !-------------------------------------------------------------------!
  ! Settings for time-integration
  !-------------------------------------------------------------------!

  ! The defaults are already set in global_varaibles.f90.
  ! The user is free to change here too.
  
!!$  dT         = 0.01_dp
!!$  start_time = 0.0_dp
!!$  end_time   = 1.0_dp
  
  aa   = 1.0_dp/dT    ! used in jacobian assembly and state update
  bb   = 1.0_dp/dt**2 ! used in jacobian assembly and state update

  time = start_time

  march: do while ( time .le. end_time)  ! loop for time-marching

     time = time + dT ! update the time

     ! Newton iterations
     newton: do k = 1, MAX_NEWTON_ITER

        ! (3) Assemble the residual and jacobian using the body
        
        ! call this method to update the state for every other iter
        if (k .gt. 1) then
           call set_state(q, q_dot, body1)
        end if

        ! residual assembly
        call disp(" >> Assembling the residuals...")
        res  = get_residual(body1)
        call disp(" >> Residual assembly complete...")
        call disp("   R   =   ", res)

        ! jacobian assembly
        call disp(" >> Assembling the Jacobian...")
        jac  = get_jacobian(body1)
        call disp(" >> Jacobian assembly complete...")
        call disp("   JAC =   ",jac)

        ! ******************************************************
        ! (3) Solve the linear system and compute update delta_q
        ! ******************************************************

        call disp(" >> Calling the linear solver...") 

        !dq = linear_solve(jac,-res,'GMRES')

        call disp(" >> Update to residual is computed...") 

        l2_norm   = norm2(dq); linf_norm =  maxval(abs(dq));

        call disp(" >> L2-Norm of the update    :", l2_norm) 
        call disp(" >> Infty-Norm of the update :", linf_norm)

        !-------------------------------------------------------------!
        ! (4) Update the state variables and then update the body
        !-------------------------------------------------------------!

        q            = get_updated_q(q, dq)
        q_dot        = get_updated_q_dot(q_dot, dq)
        !  q_double_dot = get_updated_q_double_dot(q_double_dot, dq)

        call disp("   q      =   ", q)
!        call disp("   qdot   =   ", q_dot)

        if (l2_norm .le. REL_TOL .OR. linf_norm.le. ABS_TOL) &
             & exit newton

        if (k .eq. MAX_NEWTON_ITER) then
           call disp(" >> Newton solution failed in ", k , "iterations")
           stop
        end if

     end do newton

     ! print the summary

  end do march

  call disp(" >> End of program...") 
  call disp("========================================================")
  
end program pendulum

!*********************************************************************|
! Example implementation of the function call to compute FD derivative
!*********************************************************************|
subroutine residual(q, f)

  use global_constants, only: sp, dp, TOT_NDOF

  real(dp), intent(in), dimension(TOT_NDOF)  :: q
  real(dp), intent(out), dimension(TOT_NDOF) :: f

  !  f = get_residual();  

  f = (/ sin(q(1))*sin(q(3)), sin(q(2)), sin(q(3)), sin(q(4)),&
       & sin(q(5)), sin(q(6)), &
       &sin(q(7)), sin(q(8)), sin(q(9)), &
       &sin(q(10)), sin(q(11)), sin(q(12)*sin(q(1))) /)

end subroutine residual

