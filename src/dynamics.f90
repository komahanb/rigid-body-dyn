!=====================================================================!
! Module to that contains routines to simulate the dynamics of the 
! system
!=====================================================================!

module dynamics

  use global_variables, only: solver_type, ierr
  use dynamics_input, only: read_input
  use global_constants, only: dp
  use types
  implicit none

contains

  subroutine setup_dynamics()

    use dynamics_input

    ! read the input file
    call read_input()

    ! any post processing may be?

  end subroutine setup_dynamics

  !*******************************************************************!
  ! Routine that wraps the whole execution logic from creating geometry
  ! to solving the system to time-marching
  !*******************************************************************!

  subroutine execute()

    !-----------------------------------------------------------------!
    ! setup the problem
    !-----------------------------------------------------------------!

    call setup_dynamics()

    !-----------------------------------------------------------------!
    ! solve the problem
    !-----------------------------------------------------------------!

    select case(trim(solver_type))

    case('newton-krylov')

       ! call newton_krylov_solve()

    case('LU')

       ! call lower_upper_solve()

    case DEFAULT

       !  call direct_solve()

    end select


  end subroutine execute

end module dynamics
