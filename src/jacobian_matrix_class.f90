module jacobian_matrix_class

  ! module references
  use sparse_matrix_class,  only: sparse_matrix

  ! module options
  implicit none
  private
  public :: jacobian_matrix

  ! module definitions
  type, extends(sparse_matrix) :: jacobian_matrix

     ! type options
  private

  ! type attributes
  ! pointer to a material array object
  ! pointer to a geometry object 

end type jacobian_matrix

! interfaces
interface jacobian_matrix
  procedure constructor ! add constructor to petsc_matrix generic interface
end interface jacobian_matrix

! procedures
contains

  !===================================================================!
  ! initializes an instance of jacobian_matrix
  !===================================================================!

function constructor(row_size, col_size) result(this)

  ! arguments
 integer :: col_size
 integer :: row_size
 type(jacobian_matrix) :: this

 ! begin execution

 ! set up size of matrix
 call this % set_row_size(row_size)
 call this % set_col_size(col_size) 

 ! set up preallocation

 ! associate material and geometry objects

 print *, 'I AM CONSTRUCTING'

end function constructor

end module jacobian_matrix_class
