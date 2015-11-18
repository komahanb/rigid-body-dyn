!=====================================================================!
! Module that solves the system using newton_iterations for the next
! update
!=====================================================================!

module system_solve
  
  use global_constants, only : dp
  
  use global_variables, only: solver_type, filenum, filename, &
       & ABS_TOL, REL_TOL, MAX_NEWTON_ITER, &
       & start_time, dT, time, end_time, &
       & aa, bb, fcnt, unsteady,&
       & update_norm, res_norm, newton_cnt, master,&
       & dynsys, newton_success

  use dispmodule, only: disp

  implicit none

  private
  public  :: steady_solve, unsteady_solve

contains

  !*******************************************************************!
  ! Routine that solves the system to steady state starting from the 
  ! initial guess of state values provided
  !*******************************************************************!

  subroutine steady_solve

    ! local variables
 
    integer :: k  ! loop variable
    
       !--------------------------------------------------------------!
       ! Newton iteration loop
       !--------------------------------------------------------------!

       newton_cnt = 0

       newton: do k = 1, MAX_NEWTON_ITER

          newton_cnt = newton_cnt + 1 ! increment the iter cnt

          newton_success = .false.  ! reset before starting a new iter

          !-----------------------------------------------------------!
          ! Update the state for every other iteration that first
          !-----------------------------------------------------------!

          if (k .gt. 1)  call set_state(q, q_dot, dynsys%bodies(1))
          
!!$          !-------------------------------------------------------------!
!!$          !---------------------RESIDUAL ASSEMBLY-----------------------!
!!$          !-------------------------------------------------------------!
!!$          res  = get_residual(body1)
!!$
!!$          !-------------------------------------------------------------!
!!$          !---------------------JACOBIAN ASSEMBLY-----------------------!
!!$          ! (A) Actual jacobian
!!$          ! (B) Finite difference approxiamtion to Jacobian
!!$          !-------------------------------------------------------------!
!!$
!!$          !jac = get_jacobian(body1)  ! actual jacobian
!!$
!!$          jac = finite_difference2(q, q_dot, aa, 1.0d-6) !finite diff
!!$
!!$          !-------------------------------------------------------------!
!!$          !--------------------SOLUTION TO LINEAR SYSTEM ---------------!
!!$          ! (a) Direct solution
!!$          ! (b) LU decomposition
!!$          ! (c) Iterative solution (should implement)
!!$          !-------------------------------------------------------------!
!!$
!!$          jac = direct_solve(jac); dq = matmul(jac, -res);
!!$
!!$          !dq = direct_solve(jac, -res, TOT_NDOF)
!!$
!!$          !-------------------------------------------------------------!
!!$          ! Calcualte residual tolratances                    
!!$          !-------------------------------------------------------------!
!!$
!!$          update_norm  =  maxval(abs(dq)); res_norm  =  maxval(abs(res));
!!$
!!$          if (k .eq. 1)   write(filenum,'(i4, 26F25.16)') &
!!$               & newton_cnt, &
!!$               & update_norm, res_norm, &
!!$               & (q(i), i = 1, size(q)), &
!!$               & (q_dot(i),i = 1, size(q_dot))
!!$
!!$          !-------------------------------------------------------------!
!!$          ! Update the state variables and then update the body
!!$          !-------------------------------------------------------------!
!!$
!!$          q            = get_updated_q(q, dq)
!!$          q_dot        = get_updated_q_dot(q_dot, dq)
!!$          !q_double_dot = get_updated_q_double_dot(q_double_dot, dq) 
!!$
!!$          !-------------------------------------------------------------!
!!$          ! Check for convergence and stop is necessary
!!$          !-------------------------------------------------------------!
!!$
!!$          !is converged?
!!$          if ( update_norm .le. ABS_TOL .AND. res_norm .le. ABS_TOL) then
!!$
!!$             exit newton
!!$          end if
!!$
!!$          ! is disverged?
!!$          if ( update_norm .ge. 1.0d5 .AND. res_norm .le. 1.0d5) then           
!!$             call disp(" > Solution diverging - aborting time integrtn")
!!$             exit newton
!!$          end if
!!$
!!$          !? max iters reached
!!$          if (k .eq. MAX_NEWTON_ITER) then
!!$             call disp(" >> Newton solution failed in ", k , "iterations")
!!$             call print_body(body1)
!!$             stop
!!$          end if
!!$
       end do newton


     end subroutine steady_solve

  
  !*******************************************************************!
  ! Routine that performs time integration of the system starting 
  ! from intial state to a final state of time
  !*******************************************************************!
  
  subroutine unsteady_solve

    !-----------------------------------------------------------------!
    ! Time integration loop
    !-----------------------------------------------------------------!

    write(*, *)'          Time' , '       ||dq||', '       ||R||',&
         &'       Niter',' FCNT'
    
    time = start_time
    
    time_march: do while ( time .le. end_time)
       
       ! Solve to steady state at the current time
       call steady_solve()

       ! Print the summary of the newton iteration
       write(*,'(f15.2, e15.6, e15.6, i4, i4)') &
            & time, update_norm, res_norm, newton_cnt, fcnt
       
       ! Update the time
       time = time + dT
       
    end do time_march
    
  end subroutine unsteady_solve

end module system_solve
