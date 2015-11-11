!===========================================================
! program to run dynamics simulation of a physical pendulum
!===========================================================
program pendulum

  ! import the necessary modules
  
  use global_constants, only: sp, dp, NUM_SPAT_DIM, TOT_NDOF
  use global_variables,only: ABS_TOL, REL_TOL, start_time, dT, &
       &end_time, aa, bb, MAX_NEWTON_ITER, fcnt, init
  use utils,only:deg2rad
  use types, only: body, vector
  use solver_utils, only: get_updated_q, get_updated_q_dot,&
       & get_updated_q_double_dot
  use bodies, only: create_body, set_state
  use residual, only: get_residual
  use jacobian, only: get_jacobian
  use finite_diff, only:finite_difference2
  use linear_system,only:direct_solve
  use bodies,only:print_body, create_body
  use dispmodule,only:disp
 
  implicit none

  ! position of the point of interest on the body measured in body frame
  real(dp), dimension(NUM_SPAT_DIM)  :: re    

  ! mass of the body
  real(dp)                           :: mass 

  ! the state
  real(dp)                           :: q(TOT_NDOF)= 0.0_dp

  ! time derivative of state
  real(dp)                           :: q_dot(TOT_NDOF) = 0.0_dp

  ! residual
  real(dp)                           :: res(TOT_NDOF)

  ! jacobian
  real(dp)                           :: jac(TOT_NDOF, TOT_NDOF)

  ! update from newton iterations
  real(dp)                           :: dq(TOT_NDOF) = 0.0_dp

  ! body -- single pendulum
  type(body)                         :: body1

  ! time, and norms to stop integration
  real(dp)                           :: time, update_norm, res_norm

  ! loop variable
  integer(sp)                        :: k, newton_cnt

  call disp("========================================================")
  call disp("'                 Rigid body dynamics                  '")
  call disp("========================================================")

  !-------------------------------------------------------------------!
  ! Init routine that has sanity check on the input settings
  !-------------------------------------------------------------------!
  ! Note: 
  !
  ! Has nothing as of now, we can use to:
  ! Load input files, external calls etc
  !-------------------------------------------------------------------!
  call init()

  !-------------------------------------------------------------------!
  ! set the initial states and attributes of the body
  !-------------------------------------------------------------------!
  call disp(" >> Setting initial state of the body...")

  !call random_seed(); call random_number(q);call random_number(q_dot);

  mass = 2.0_dp
  re   = (/ 0.1_dp, 0.2_dp, 0.3_dp /)

  ! define the initial states
  q(1:3)    = (/ 1.0_dp, 2.0_dp, 3.0_dp /)
  q(4:6)    = (/ deg2rad(20.0_dp), deg2rad(0.0_dp), deg2rad(20.0_dp)/)
  q(7:9)    = (/ 1.0_dp,  2.0_dp, 1.0_dp /)
  q(10:12)  = (/ 2.0_dp, 4.0_dp, 5.0_dp /)

  ! define the intial time derivative of state
  q_dot(1:3)   = (/ 2.0_dp, 2.0_dp, 2.0_dp /)
  q_dot(4:6)   = (/ 0.25_dp, 0.0_dp, 0.0_dp /)
  q_dot(7:9)   = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
  q_dot(10:12) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
  !-------------------------------------------------------------------!
  ! create a pendulum body using the state and attr above
  !-------------------------------------------------------------------!

  call disp(" >> Creating a body...")

  body1 = create_body(mass, vector(re), q, q_dot)

  call print_body(body1)


  !  call disp("   q0  =   ", q)
  !  call disp("   qdot0  =   ", q_dot)

  !-------------------------------------------------------------------!
  ! Settings for time-integration
  !-------------------------------------------------------------------!

!!$! The following defaults are already set in global_varaibles.f90.
!!$! The user is free to change here too.
!!$  dT         = 0.01_dp
!!$  start_time = 0.0_dp
!!$  end_time   = 1.0_dp

  aa   = 1.0_dp/dT    ! used in jacobian assembly and state update
  bb   = 1.0_dp/dT**2 ! used in jacobian assembly and state update
  time = start_time

  ! print an initial summary
  call disp(">> Starting time-integration...")
  call disp("  >Initial time        :", start_time)
  call disp("  >End time            :", end_time)
  call disp("  >Time-step           :", dT)
  call disp("  >Jacobian co-eff     :", aa)

  call disp("  >Max Newton iters    :", MAX_NEWTON_ITER)
  call disp("  >Relative tolenrance :", REL_TOL)
  call disp("  >Absolute tolenrance :", ABS_TOL)

  print*, '          time' , '       ||dq||', '       ||R||',&
       &  '             KE', '          PE','               TE', &
       &'       Niter',' FCNT'

  march: do while ( time .le. end_time)  ! loop for time-marching

     time = time + dT ! update the time

     ! Newton iterations
     newton_cnt = 0

     newton: do k = 1, MAX_NEWTON_ITER

        newton_cnt = newton_cnt + 1

        ! call this method to update the state for every other iter
        ! we don't have to create a whole new body again
        if (k .gt. 1) then
           call set_state(q, q_dot, body1)
           !stop
        end if

        !call disp(" >> Assembling the residuals...")
        res  = get_residual(body1)
        !call disp(" >> Residual assembly complete...")
        !jac = get_jacobian(body1) 
        jac = finite_difference2(q, q_dot, aa, 1.0d-6) !finite diff

        !-------------------------------------------------------------!
        ! Solve the linear system and compute update delta_q
        !-------------------------------------------------------------!
        !call disp(" >> Calling the linear solver...") 
        !dq = linear_solve(jac,-res,'GMRES')
        dq = direct_solve(jac, -res, TOT_NDOF)
        !call disp(" >> Update to state residual is computed...") 

        update_norm   =  maxval(abs(dq));
        res_norm =  maxval(abs(res));

        !       call disp(" >> L2-Norm of the update    :", update_norm) 
        !       call disp(" >> Infty-Norm of the update :", res_norm)

        !-------------------------------------------------------------!
        ! Update the state variables and then update the body
        !-------------------------------------------------------------!

        q            = get_updated_q(q, dq)
        q_dot        = get_updated_q_dot(q_dot, dq)
        ! q_double_dot = get_updated_q_double_dot(q_double_dot, dq) 

        !-------------------------------------------------------------!
        ! Check for convergence and stop is necessary
        !-------------------------------------------------------------!
        if ( update_norm .le. REL_TOL .AND. res_norm .le. REL_TOL) &
             & exit newton

        if (k .eq. MAX_NEWTON_ITER) then
           call disp(" >> Newton solution failed in ", k , "iterations")
           stop
        end if

     end do newton

     ! print the summary of the newton iteration
     write(*,'(f15.2,e15.6,e15.6,e15.6,e15.6,e15.6,xi4,xi4)') &
          &time, update_norm, res_norm, body1%KE, body1%PE, &
          & body1%KE + body1%PE, newton_cnt, fcnt

  end do march

  call print_body(body1)

  call disp(" >> End of program...") 

end program pendulum

!*********************************************************************!
! Callbaclk function for approximation of jacobian using finite-diff
!*********************************************************************!
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

  ! we use these params in this problem
  mass  = 2.0_dp;   re   = (/ 0.1_dp, 0.2_dp, 0.3_dp /);

  body1 = create_body(mass, vector(re), q, qdot);

  f = get_residual(body1);  


end subroutine residual2

!*********************************************************************!
! Callbaclk function for approximation of jacobian using finite-diff
!*********************************************************************!
subroutine residual(q, f)

  use global_constants, only: sp, dp, TOT_NDOF

  real(dp), intent(in), dimension(TOT_NDOF)  :: q
  real(dp), intent(out), dimension(TOT_NDOF) :: f

  f = (/ sin(q(1))*sin(q(3)), sin(q(2)), sin(q(3)), sin(q(4)),&
       & sin(q(5)), sin(q(6)), &
       &sin(q(7)), sin(q(8)), sin(q(9)), &
       &sin(q(10)), sin(q(11)), sin(q(12)*sin(q(1))) /)

  stop"dummy impl"

end subroutine residual

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
