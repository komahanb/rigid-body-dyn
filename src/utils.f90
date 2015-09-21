module utils
use constants

implicit none

! overload * for cross product
interface operator (*)
   module procedure cross_pdt, scal_vec, vector_matrix, matrix_vector, matrix_matrix
end interface operator (*)

! overload abs 
interface abs
   module procedure abs_vec
end interface abs

! constructor for vector data type
interface vector
   module procedure new_vec, get_array
end interface vector


! constructor for vector data type
interface matrix
   module procedure new_matrix, get_matrix
end interface 

contains


function scal_vec(a, b) 
  real(dp), intent (in)      :: a
  type (vector), intent (in) :: b
  type (vector) :: scal_vec
  scal_vec%x = a*b%x
  scal_vec%y = a*b%y
  scal_vec%z = a*b%z
end function scal_vec


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

  matrix_matrix = new_matrix_from_matrix( matmul( get_matrix(A), get_matrix(B) ) )

end function matrix_matrix





! returns a skew-symmetric matrix for doing cross product
function skew(a) result(tilde_a)
  real(dp), intent(in)           :: a(3)
  real(dp)                       :: tilde_a(3,3)
  tilde_a=reshape((/ 0.0_dp, a(3), -a(2), -a(3), 0.0_dp, a(1),  a(2), -a(1), 0.0_dp /), (/3, 3/))
end function skew

! returns cross product of vectors a and b
function cross(a ,b) 
  real(dp), intent(in) :: a(3), b(3)
  real(dp)             :: cross(3)
  cross =matmul(skew(a), b) 
end function cross

! a different implementation of cross product
function cross_pdt(a, b) 
  type (vector), intent (in) :: a, b
  type (vector) :: cross_pdt
  cross_pdt%x = a%y * b%z - a%z * b%y
  cross_pdt%y = a%z * b%x - a%x * b%z
  cross_pdt%z = a%x * b%y - a%y * b%x
end function cross_pdt

! returns the dot product of two vectors
function dot(a,b)
  type (vector), intent (in) :: a, b
  real(dp)                   :: dot
  dot = a%x*b%x + a%y*b%y + a%z*b%z
end function dot

! returns the magnitude of a vector
function abs_vec(a)
  type (vector), intent (in) :: a
  real(dp)                   :: abs_vec
    abs_vec = sqrt(a%x**2 + a%y**2 + a%z**2)
end function abs_vec

! constructor for a new vector
function new_vec(a)
  real(dp), intent(in) :: a(3)
  type(vector) :: new_vec
  
  new_vec%x=a(1)
  new_vec%y=a(2)
  new_vec%z=a(3)

end function new_vec

! get the vector entries as array
function get_array(a)
  type(vector), intent(in) :: a
  real(dp)                 :: get_array(num_spat_dim)

  get_array(1) = a%x
  get_array(2) = a%y
  get_array(3) = a%z

end function get_array

! constructor for a new matrix
function new_matrix(a)
  real(dp), intent(in) :: a(num_spat_dim**2)
  type(matrix) :: new_matrix
  real(dp)     :: tmp(num_spat_dim, num_spat_dim)

  tmp = reshape(a, (/ num_spat_dim, num_spat_dim /))

  new_matrix%ij = tmp

end function new_matrix

! constructor for a new matrix
function new_matrix_from_matrix(A)
  real(dp), intent(in) :: A(num_spat_dim, num_spat_dim)
  type(matrix)         :: new_matrix_from_matrix
  new_matrix_from_matrix%ij=A
end function new_matrix_from_matrix

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



end module utils
