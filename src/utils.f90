module utils
use constants

implicit none

! overload * for cross product
interface operator (*)
   module procedure cross, scal_vec, scal_matrix, vector_matrix, matrix_vector, matrix_matrix
end interface operator (*)

interface operator (+)
   module procedure add_matrices
end interface operator (+)

interface operator (-)
   module procedure sub_matrices, negate_vector,  negate_matrix
end interface operator (-)




! overload abs 
interface abs
   module procedure abs_vec
end interface abs

! constructor for vector data type
interface vector
   module procedure get_array
end interface vector

! constructor for vector data type
interface matrix
   module procedure new_matrix_from_array, get_matrix
end interface 

contains

function scal_vec(a, b) 
  real(dp), intent (in)      :: a
  type (vector), intent (in) :: b
  type (vector) :: scal_vec
  scal_vec%x   = a*b%x
end function scal_vec

! {aB} = a[B]
function scal_matrix(a, B) 
 
  real(dp), intent (in)     :: a
  type(matrix), intent (in) :: B 
  type(matrix)              :: scal_matrix

  scal_matrix%ij =  a*B%ij

end function scal_matrix

! {c} = {a}^T[B]
function vector_matrix(a, B) 
  
  type(vector), intent (in) :: a
  type(matrix), intent (in) :: B 
  type(vector) ::  vector_matrix

end function vector_matrix

! {c} = [B]{a}
function matrix_vector(B, a) 
  
  type(vector), intent (in) :: a
  type(matrix), intent (in) :: B 
  type(vector) ::  matrix_vector

end function matrix_vector

! [C] = [A]{B}
function matrix_matrix(A, B) 

  type(matrix), intent (in) :: A, B 
  type(matrix)              :: matrix_matrix

  matrix_matrix = matrix( matmul( get_matrix(A), get_matrix(B) ) )

end function matrix_matrix

! returns a skew-symmetric matrix for doing cross product
function skew(a) result(tilde_a)
  type(vector), intent(in)           :: a
  type(matrix)                       :: tilde_a
  tilde_a = matrix((/ 0.0_dp, a%x(3), -a%x(2), -a%x(3), 0.0_dp, a%x(1),  a%x(2), -a%x(1), 0.0_dp /))
end function skew

! returns cross product of vectors a and b
function cross(a ,b) 
  type(vector), intent(in) :: a, b
  type(vector)             :: cross
!  real(dp) ::Amat(3,3)
!  real(dp) :: x(3)
!  Amat = get_matrix(skew(a))
!  x = get_array(b)

  cross = vector(matmul(get_matrix(skew(a)), get_array(b)))
 
end function cross

! returns the dot product of two vectors
function dot(a,b)
  type (vector), intent (in) :: a, b
  real(dp)                   :: dot
  dot = norm2(a%x*b%x)
end function dot

! returns the dot product of two vectors
function square(a)
  type(vector), intent (in) :: a
  type(vector)              :: square
  square = vector(a%x(:)**2)
end function square


! returns the magnitude of a vector
function abs_vec(a)
  type (vector), intent (in) :: a
  real(dp)                   :: abs_vec  
  abs_vec = norm2(a%x(:)**2)
end function abs_vec

! constructor for a new vector
!!$function new_vec(a)
!!$  real(dp), intent(in) :: a(num_spat_dim)
!!$  type(vector) :: new_vec
!!$  new_vec%x=a
!!$end function new_vec

! get the vector entries as array
function get_array(a)
  type(vector), intent(in) :: a
  real(dp)                 :: get_array(num_spat_dim)
  get_array=a%x
end function get_array

! constructor for a new matrix
function new_matrix_from_array(a)
  real(dp), intent(in) :: a(num_spat_dim**2)
  type(matrix) :: new_matrix_from_array
  new_matrix_from_array%ij = reshape(a, (/ num_spat_dim, num_spat_dim /))
end function new_matrix_from_array

! constructor for a new matrix
function mat(arr,n)
  real(dp), intent(in)    :: arr(n*n)
  integer(sp), intent(in) :: n
  real(dp)                :: mat(n,n)
  mat   = reshape(arr, (/n,n/))
end function mat

!!$! constructor for a new matrix
!!$function new_matrix_from_matrix(A)
!!$  real(dp), intent(in) :: A(num_spat_dim, num_spat_dim)
!!$  type(matrix)         :: new_matrix_from_matrix
!!$  new_matrix_from_matrix%ij=A
!!$end function new_matrix_from_matrix

!!$! constructor for a new matrix
!!$function new_matrix(a,m,n)
!!$  integer(sp), intent(in) :: m, n
!!$  real(dp), intent(in)    :: a(m*n)
!!$  type(matrix)            :: new_matrix
!!$  real(dp)                :: tmp(num_spat_dim, num_spat_dim)
!!$
!!$  tmp = reshape(a, (num_spat_dim, num_spat_dim /))
!!$  
!!$  new_matrix%ij = tmp
!!$
!!$end function new_matrix


! get the matrix entries as array
function get_matrix(A)
  type(matrix), intent(in) :: A
  real(dp)                 :: get_matrix(num_spat_dim, num_spat_dim)
  get_matrix = A%ij
end function get_matrix

! transpose of a matrix
function trans(A)
  type(matrix),intent(in) :: A
  type(matrix)            :: trans
  trans = matrix(transpose(get_matrix(A)))
end function trans

!returns the matrix addition of matrices of TYPE matrix
function add_matrices(A, B)
  type(matrix), intent(in) :: A, B
  type(matrix)  :: add_matrices
  add_matrices%ij = A%ij +B%ij
end function add_matrices

!returns the matrix subtraction of matrices of TYPE matrix
function sub_matrices(A, B)
  type(matrix), intent(in)  :: A, B
  type(matrix)  :: sub_matrices
  sub_matrices%ij = A%ij -B%ij
end function sub_matrices

! returns the negative of type matrix
function negate_matrix(A)
  type(matrix), intent(in)  :: A
  type(matrix)  :: negate_matrix
  negate_matrix%ij = -A%ij
end function negate_matrix

! returns the negative of type vector
function negate_vector(a)
  type(vector), intent(in)  :: a
  type(vector)  :: negate_vector
  negate_vector%x = -a%x
end function negate_vector

end module utils
