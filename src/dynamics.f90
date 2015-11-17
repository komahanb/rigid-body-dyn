!=====================================================================!
! Module to that contains routines to simulate the dynamics of the 
! system
!=====================================================================!

module dynamics

  use filehandler, only: newunit

  use global_constants, only: sp, dp

  use global_variables, only: solver_type, filenum, filename, &
       & ABS_TOL, REL_TOL, MAX_NEWTON_ITER, &
       & start_time, dT, time, end_time, &
       & aa, bb, fcnt, unsteady,&
       & update_norm, res_norm, newton_cnt, master

  use system_input_handler, only: read_system_input

  use dispmodule, only:disp

  implicit none

contains

  !*******************************************************************!
  ! Set up the system including
  ! (a) reading the input
  ! (b) any pre-processing
  ! (c) 
  !*******************************************************************!

  subroutine setup_system()

    ! read the input file
    call read_system_input()

    ! any pre-processing may be necessary?

    call read_solver_input()

    ! open solution.inp
    ! call setup_integration() dt, time_step, start_time, end_time

  end subroutine setup_system

  !*******************************************************************!
  ! Reads the input file for solution method and sets up  global vars
  !*******************************************************************!

  subroutine read_solver_input()

    call disp(" >> Loading solver settings from solution.inp") 
    
    ! get a new unit number
    filenum = newunit()

    ! read the file
    open (unit=filenum, file="solution.inp")
    
    read(filenum, *) start_time
    read(filenum, *) end_time
    read(filenum, *) dT
    
    read(filenum, *) MAX_NEWTON_ITER
    read(filenum, *) REL_TOL
    read(filenum, *) ABS_TOL
    
    close(filenum)
    
    ! steady or unsteady?
    if (start_time .ne. end_time) unsteady = .true.

    ! a few sanity checks for unsteady 
    if (unsteady) then

       if (start_time .gt. end_time)&
            &stop "ERROR: Start time is greater than end time"

       if (dT .lt. epsilon(1.0_dp)) &
            &stop "ERROR: Time step less than machine precision"

    end if

    ! everything fine print summary
    call disp("")
    call disp(" >> Solver inputs are loaded successfully ") 
    call disp("  >Initial time        :", start_time)
    call disp("  >End time            :", end_time)
    call disp("  >Time-step           :", dT)
    call disp("  >Max Newton iters    :", MAX_NEWTON_ITER)
    call disp("  >Relative tolenrance :", REL_TOL)
    call disp("  >Absolute tolenrance :", ABS_TOL)
    call disp("  >Unsteady problem    :", unsteady)
    call disp("")

  end subroutine read_solver_input

  !*******************************************************************!
  ! solve the system using appropriate linear solver
  !*******************************************************************!

  subroutine solve_system()

    !-------------------------------------------------------------------!
    ! Settings for time-integration
    !-------------------------------------------------------------------!

    aa   = 1.0_dp/dT    ! used in jacobian assembly and state update
    bb   = 1.0_dp/dT**2 ! used in jacobian assembly and state update

    call disp("  >Jacobian co-eff alpha:", aa)
    call disp("  >Jacobian co-eff beta:" , bb)

    call disp("")
    call disp(">> Starting time-integration...")
    
    print*, '          time' , '       ||dq||', '       ||R||',&
         !         &  '             KE', '          PE','               TE', &
         &'       Niter',' FCNT'

    !-----------------------------------------------------------------!
    ! Open up a file unit to record the results       
    !-----------------------------------------------------------------!

    filenum = newunit()

    open(unit=filenum,file="state.dat",action="write",status="replace")

    ! Write tecplot labels for state vars and their time derivs
    write(filenum,*) 'VARIABLES = "ITER" "R1" "R2" "R3" "T1" "T2" &
         &"T3" "V1" "V2" "V3" "W1" "W2" "W3" "DR1" "DR2" "DR3" &
         & "DT1" "DT2" "DT3" "DV1" "DV2" "DV3" "DW1" "DW2" "DW3" "NORM"'

    !-------------------------------------------------------------------!
    ! Time integration loop
    !-------------------------------------------------------------------!

    time = start_time

    time_march: do while ( time .le. end_time)  ! loop for time-marching

!!$       !----------------------------------------------------------------!
!!$       ! Newton iteration loop
!!$       !----------------------------------------------------------------!
!!$
!!$       newton_cnt = 0
!!$
!!$       newton: do k = 1, MAX_NEWTON_ITER
!!$
!!$          newton_cnt = newton_cnt + 1 ! increment the iter cnt
!!$
!!$          newton_success = .false.  ! reset before starting a new iter
!!$
!!$          !-------------------------------------------------------------!
!!$          ! Update the state for every other iteration that first
!!$          !-------------------------------------------------------------!
!!$          if (k .gt. 1)  call set_state(q, q_dot, body1)
!!$
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
!!$       end do newton

       !--------------------------------------------------------------!
       ! print the summary of the newton iteration
       !--------------------------------------------------------------!
       write(*,'(f15.2, e15.6, e15.6, xi4, xi4)') &
            & time, update_norm, res_norm, newton_cnt, fcnt

       time = time + dT ! update the time

    end do time_march
    
    ! close the file that records the state vars across time
    
    close (filenum)

!!$    select case(trim(solver_type))
!!$
!!$    case('newton-krylov')
!!$
!!$       ! call newton_krylov_solve()
!!$
!!$    case('LU')
!!$
!!$       ! call lower_upper_solve()
!!$
!!$    case DEFAULT
!!$
!!$       !  call direct_solve()
!!$
!!$    end select

  end subroutine solve_system


  !*******************************************************************!
  ! Post process the results of the simulations 
  ! e.g. tecplot output
  !      iteration history
  !*******************************************************************!

  subroutine post_process()

    ! print the system here at its final state
    !    call print_body(body1)

    call disp(" >> End of program...") 

  end subroutine post_process

  !*******************************************************************!
  ! Routine that wraps the whole execution logic from creating geometry
  ! to solving the system to time-marching
  !*******************************************************************!

  subroutine execute()

    !-----------------------------------------------------------------!
    !------------------- SETUP THE PROBLEM----------------------------!
    !-----------------------------------------------------------------!

    call setup_system()

    !-----------------------------------------------------------------!
    !------------------- SOLVE THE SYSTEM ----------------------------!
    !-----------------------------------------------------------------!


    call solve_system()

    !-----------------------------------------------------------------!
    !-------------------- POST-PROCESS THE RESULTS -------------------!
    !-----------------------------------------------------------------!

    call post_process()

  end subroutine execute

end module dynamics
