!=====================================================================!
! Interface to Newton-type, quasi-Newton, full approximation scheme 
! (FAS) multigrid, and other methods for solving systems of nonlinear 
! equations in PETSc
!=====================================================================!
module snes

  ! module options
  implicit none
  private
  public solver

  interface solver
     module procedure solve
  end interface solver
  
contains

  subroutine solve


  end subroutine solve

end module snes
