module utils
use constants

implicit none


! overload * for cross product
interface operator (*)
   module procedure cross, scal_vec, scal_matrix, vector_matrix, matrix_vector, matrix_matrix
end interface operator (*)

interface operator (+)
   module procedure add_matrices, add_vectors
end interface operator (+)

interface operator (-)
   module procedure sub_matrices, sub_vectors, negate_vector, negate_matrix
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

type(vector), parameter       :: zeroV = vector((/ dzero, dzero, dzero /))
type(vector), parameter       :: unitV = vector((/ dzero, dzero, dzero /))

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

! get the vector entries as array
function get_array_1d(a,n)
  integer :: n
  type(vector), intent(in) :: a(n)
  real(dp)                 :: get_array_1d(n*num_spat_dim)
  integer                  :: j, is_j, ie_j
  do j = 1, n
     call split(j,is_j,ie_j) ! split j index storage
     get_array_1d(is_j:ie_j) = get_array(a(j))
  end do
end function get_array_1d

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

! unwraps a 2d matrix and stores as real numbers
function get_matrix_2d(A,m,n)
  integer(sp) :: m,n
  type(matrix), intent(in) :: A(m,n)
  real(dp)    :: get_matrix_2d(m*num_spat_dim, n*num_spat_dim) 
  integer(sp) :: i,j
  integer(sp) :: is_i, ie_i, is_j, ie_j

  do i = 1, n
     do j = 1, m
        call split(j,is_j,ie_j) ! split j index storage
        call split(i,is_i,ie_i) ! split i index storage
        get_matrix_2d(is_j:ie_j,is_i:ie_i) = get_matrix(A(j,i))
        !        print*, is_j,ie_j,is_i,ie_i,i,j
     end do
  end do

end function get_matrix_2d

 subroutine split(i,is,ie)
   integer(sp) :: i
   integer(sp ):: is, ie
   is = (num_spat_dim*(i-1) ) + 1
   ie = num_spat_dim*( (i-1) + 1)
!   print*, is,ie !,num_spat_dim,i
 end subroutine split

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



!returns the matrix addition of vectors of TYPE matrix
function add_vectors(a, b)
  type(vector), intent(in) :: a, b
  type(vector)  :: add_vectors
  add_vectors%x = a%x + b%x
end function add_vectors


!returns the matrix subtraction of matrices of TYPE matrix
function sub_matrices(A, B)
  type(matrix), intent(in)  :: A, B
  type(matrix)  :: sub_matrices
  sub_matrices%ij = A%ij -B%ij
end function sub_matrices



!returns the matrix addition of vectors of TYPE matrix
function sub_vectors(a, b)
  type(vector), intent(in) :: a, b
  type(vector)  :: sub_vectors
  sub_vectors%x = a%x - b%x
end function sub_vectors


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

! convert from degree to radian
function deg2rad(deg)
  real(dp) :: deg2rad,deg
  deg2rad  =  deg/rad_to_deg
end function deg2rad

! convert from radian to degree
function rad2deg(rad)
  real(dp) :: rad2deg,rad
  rad2deg  =  rad*rad_to_deg
end function rad2deg

end module utils
