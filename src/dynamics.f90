!=====================================================================!
! Module to that contains routines to simulate the dynamics of the 
! system
!=====================================================================!

module dynamics

  use global_variables, only: solver_type
  use system_input_handler, only: setup_system

  implicit none

contains

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


  !*******************************************************************!
  ! Post process the results of the simulations 
  ! e.g. tecplot output
  !      iteration history
  !*******************************************************************!
  
  subroutine post_process()
    
  end subroutine post_process


  !*******************************************************************!
  ! solve the system using appropriate linear solver
  !*******************************************************************!
  
  subroutine solve_system()

    select case(trim(solver_type))

    case('newton-krylov')

       ! call newton_krylov_solve()

    case('LU')

       ! call lower_upper_solve()

    case DEFAULT

       !  call direct_solve()

    end select

  end subroutine solve_system

end module dynamics
