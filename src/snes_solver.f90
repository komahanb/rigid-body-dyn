module cmfd_snes_solver

  use cmfd_loss_operator,     only: loss_operator,init_M_operator,             &
       &                                  build_loss_matrix,destroy_M_operator
  use cmfd_prod_operator,     only: prod_operator,init_F_operator,             &
       &                                  build_prod_matrix,destroy_F_operator
  use cmfd_jacobian_operator, only: jacobian_operator,init_J_operator,         &
       &                                  build_jacobian_matrix,destroy_J_operator,  &
       &                                  operators
  use cmfd_slepc_solver,      only: cmfd_slepc_execute

  implicit none

#include <finclude/petsc.h90>

  type(jacobian_operator) :: jac_prec
  type(operators) :: ctx

  Mat         :: jac         ! jacobian matrix
  ! Mat         :: jac_prec    ! preconditioned jacobian matrix
  Vec         :: resvec      ! residual vector
  Vec         :: xvec        ! results
  KSP         :: ksp         ! linear solver context
  PC          :: pc          ! preconditioner
  SNES        :: snes        ! nonlinear solver context
  integer     :: ierr        ! error flag

contains

  !===============================================================================
  ! CMFD_SNES_EXECUTE
  !===============================================================================

  subroutine cmfd_snes_execute()

    ! call slepc solver 
    call cmfd_slepc_execute()

    ! initialize data
    call init_data()

    ! initialize solver
    call init_solver()

    ! solve the system
    call SNESSolve(snes,PETSC_NULL,xvec,ierr)

    ! extracts results to cmfd object
    call extract_results()

    ! deallocate all slepc data
    call finalize()

  end subroutine cmfd_snes_execute
