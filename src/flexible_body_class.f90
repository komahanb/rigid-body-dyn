module flexible_body_class

  ! module references
  use body_class
  use global_constants
  use global_variables
  use types
  use utils
  use rigid_body_class
  use dispmodule, only : disp


  ! module options
  implicit none
  private
  public :: flexible_body, print_flexible_body

  ! module definitions
  type, extends(rigid_body) :: flexible_body

     !----------------------------------------------------------------!
     ! elatic state variables and time derivatives
     !----------------------------------------------------------------!
     type(vector) :: qs
     type(vector) :: qs_dot
     type(vector) :: qs_double_dot

     type(matrix) :: p              ! 
     type(matrix) :: h              ! 

     type(matrix) :: K              ! stiffness matrix
     type(matrix) :: M              ! mass matrix

     type(vector) :: f              ! elastic force

  end type flexible_body

  ! interfaces
  interface flexible_body

     procedure constructor ! add constructor

  end interface flexible_body

contains
  !*******************************************************************!
  ! create a body using the supplied parameters
  !-------------------------------------------------------------------!
  ! mass   : mass of the body
  ! c      : first moment of mass 
  ! J      : second moment of mass
  ! fr     : reaction force
  ! gr     : reaction torque
  ! q, qdot: state vector and time derivatives
  ! qddot  : second time derivative of the state (used in elastic only)
  !*******************************************************************!
function constructor(mass, c, J, fr, gr, q, qdot, qs, &
    & qs_dot, qs_double_dot, f, p, h, K, M) result(this)

  ! inputs
 real(dp), optional, intent(in) :: mass 
 real(dp), optional, intent(in) :: c(NUM_SPAT_DIM)
 real(dp), optional, intent(in) :: J(NUM_SPAT_DIM, NUM_SPAT_DIM)
 real(dp), optional, intent(in) :: fr(NUM_SPAT_DIM)
 real(dp), optional, intent(in) :: gr(NUM_SPAT_DIM)

 real(dp), optional, intent(in) :: q(NDOF_PBODY)
 real(dp), optional, intent(in) :: qdot(NDOF_PBODY)

 real(dp), optional, intent(in) :: qs(NUM_SPAT_DIM)
 real(dp), optional, intent(in) :: qs_dot(NUM_SPAT_DIM)
 real(dp), optional, intent(in) :: qs_double_dot(NUM_SPAT_DIM)

 real(dp), optional, intent(in) :: f(NUM_SPAT_DIM)

 real(dp), optional, intent(in) :: p(NUM_SPAT_DIM, NUM_SPAT_DIM)
 real(dp), optional, intent(in) :: h(NUM_SPAT_DIM, NUM_SPAT_DIM)
 real(dp), optional, intent(in) :: K(NUM_SPAT_DIM, NUM_SPAT_DIM)
 real(dp), optional, intent(in) :: M(NUM_SPAT_DIM, NUM_SPAT_DIM)

 ! input/output
 type(flexible_body):: this

! call this%constructor(this, mass, c, J, fr, gr, q, qdot)

 if (present(qs))  this%qs = vector(qs)
 if (present(qs_dot))  this%qs_dot = vector(qs_dot)
 if (present(qs_double_dot))  this%qs_double_dot = vector(qs_double_dot)

 if (present(f))  this%f = vector(f)

 if (present(p))  this%p = matrix(p)
 if (present(h))  this%h = matrix(h)
 if (present(K))  this%K = matrix(K)
 if (present(M))  this%M = matrix(M)

end function constructor

!*******************************************************************!
! routine that prints the state and properties of the flexible body
!******************************************************************!
subroutine print_flexible_body(this)

  use dispmodule

  class(flexible_body):: this

  print*,"Implement please!!"

end subroutine print_flexible_body

end module flexible_body_class
