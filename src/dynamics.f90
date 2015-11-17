!=====================================================================!
! Module to that contains routines to simulate the dynamics of the 
! system
!=====================================================================!

module dynamics

  use global_constants, only : dp
  
  use global_variables, only : unsteady, aa, bb, dT
  
  use input_handler, only : read_system_input, read_solver_input
  
  use system_solve, only : steady_solve, unsteady_solve

  use dispmodule, only: disp

  implicit none

  private

  public :: execute

contains

  !*******************************************************************!
  ! Set up the system including
  ! (a) reading the input
  ! (b) any pre-processing
  !*******************************************************************!

  subroutine setup_problem()

    ! read dynamic system related inputs
    call read_system_input()

    ! read solution realated inputs
    call read_solver_input()

  end subroutine setup_problem


  !*******************************************************************!
  ! solve the system using appropriate linear solver
  !*******************************************************************!

  subroutine solve_problem()

    !-------------------------------------------------------------------!
    ! Settings for time-integration
    !-------------------------------------------------------------------!

    aa   = 1.0_dp/dT    ! used in jacobian assembly and state update
    bb   = 1.0_dp/dT**2 ! used in jacobian assembly and state update

!!$    call disp("  > Jacobian co-eff alpha:", aa)
!!$    call disp("  > Jacobian co-eff beta:" , bb)
!!$
!!$    call disp("")
!!$
!!$
!!$    !         &  '             KE', '          PE','               TE', &


    !-----------------------------------------------------------------!
    ! Open up a file unit to record the results       
    !-----------------------------------------------------------------!
!!$
!!$    filenum = newunit()
!!$
!!$    open(unit=filenum,file="state.dat",action="write",status="replace")
!!$
!!$    ! Write tecplot labels for state vars and their time derivs
!!$    write(filenum,*) 'VARIABLES = "ITER" "R1" "R2" "R3" "T1" "T2" &
!!$         &"T3" "V1" "V2" "V3" "W1" "W2" "W3" "DR1" "DR2" "DR3" &
!!$         &"DT1" "DT2" "DT3" "DV1" "DV2" "DV3" "DW1" "DW2" "DW3" "NORM"'
    
    ! close the file that records the state vars across time    
!!$    close (filenum)

    if (unsteady) then

       call unsteady_solve()

    else

       call steady_solve()

    end if
    
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

  end subroutine solve_problem


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

    call setup_problem()

    !-----------------------------------------------------------------!
    !------------------- SOLVE THE SYSTEM ----------------------------!
    !-----------------------------------------------------------------!


    call solve_problem()

    !-----------------------------------------------------------------!
    !-------------------- POST-PROCESS THE RESULTS -------------------!
    !-----------------------------------------------------------------!

    call post_process()

  end subroutine execute

end module dynamics
