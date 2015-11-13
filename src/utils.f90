!=====================================================================!
! Module that contains important vector and matrix operations that are 
! used throughout the code
!=====================================================================!
! dot product, cross product, skew symmetric matrix
!
!
!=====================================================================!
module utils

  use global_constants
  use global_variables
  use types

  implicit none

  !*******************************************************************!
  ! overload * for dot product and other matrix-vector operations
  !-------------------------------------------------------------------!
  ! square(xx), dot(xx, xx), xx*xx      -- > all does DOT product
  ! cross product of two vectors a and b --> skew(a)*b
  ! Other scenarios where (*) is used are:
  ! 2.0*a, 2.0*A, a*A, A*a, A*A
  !****************************************************************
  interface operator (*)
     module procedure dot, scal_vec, scal_matrix, vector_matrix, &
          &matrix_vector, matrix_matrix
  end interface operator (*)

  !*******************************************************************!
  ! overload (+) for addition of VECTOR and MATRIX  data types
  !*******************************************************************!
  interface operator (+)
     module procedure add_matrices, add_vectors
  end interface operator (+)

  !*******************************************************************!
  ! overload (-) for negation of VECTOR and MATRIX  data types
  !*******************************************************************!
  interface operator (-)
     module procedure sub_matrices, sub_vectors, negate_vector, &
          &negate_matrix
  end interface operator (-)

  !*******************************************************************!
  ! overload to get the absolute value of the VECTOR data type
  !*******************************************************************!
  interface norm
     module procedure norm_vec
  end interface norm

  !*******************************************************************!
  ! Gets plain array from VECTOR data type
  !*******************************************************************!
  interface array
     module procedure get_array, get_array_1d
  end interface array

  !*******************************************************************!
  ! Gets a multidimensional array from MATRIX data type
  !*******************************************************************1
  interface matrix
     module procedure new_matrix_from_array, get_matrix, get_matrix_2d
  end interface matrix

contains

  !**********************************
  ! product of a scalar and vector
  !**********************************
  function scal_vec(a, b) 

    real(dp), intent (in)      :: a
    type (vector), intent (in) :: b
    type (vector)              :: scal_vec

    scal_vec%x   = a*b%x

  end function scal_vec

  !*******************************************
  ! product of a scalar and matrix{aB} = a[B]
  !******************************************
  function scal_matrix(a, B) 

    real(dp), intent (in)     :: a
    type(matrix), intent (in) :: B 
    type(matrix)              :: scal_matrix

    scal_matrix%ij =  a*B%ij

  end function scal_matrix

  !**********************************************
  ! product of a vector and matrix {c} = {a}^T[B]
  !**********************************************
  function vector_matrix(a, B) 

    type(vector), intent (in) :: a
    type(matrix), intent (in) :: B 
    type(vector) ::  vector_matrix

    vector_matrix = vector(matmul(a%x, B%ij))

  end function vector_matrix

  !**********************************************
  ! product of a matrix with vector {c} = [B]{a}
  !**********************************************
  function matrix_vector(B, a) 

    type(vector), intent (in) :: a
    type(matrix), intent (in) :: B 
    type(vector) ::  matrix_vector

    matrix_vector = vector(matmul(B%ij, a%x))

  end function matrix_vector

  !**********************************************
  ! product of two matrices [C] = [A]{B}
  !**********************************************
  function matrix_matrix(A, B) 

    type(matrix), intent (in) :: A, B 
    type(matrix)              :: matrix_matrix

    matrix_matrix = matrix( matmul(A%ij,B%ij))

  end function matrix_matrix

  !*********************************************************
  ! returns a skew-symmetric matrix for doing cross product
  !*********************************************************
  function skew(a)

    type(vector), intent(in)           :: a
    type(matrix)                       :: skew

    skew = matrix((/ 0.0_dp, a%x(3), -a%x(2), -a%x(3), 0.0_dp, a%x(1),&
         &  a%x(2), -a%x(1), 0.0_dp /))

  end function skew

  !*************************************************************
  ! returns cross product of vectors a and b (can use skew too)
  !*************************************************************
  function cross(a ,b) 

    type(vector), intent(in) :: a, b
    type(vector)             :: cross

    cross = vector(matmul(get_matrix(skew(a)), get_array(b)))

  end function cross

  !****************************************
  ! returns the dot product of two vectors
  !****************************************
  function dot(a,b)

    type (vector), intent (in) :: a, b
    real(dp)                   :: dot

    dot = sum(a%x*b%x)

  end function dot

  !*****************************************************
  ! returns the magnitude of a vector
  !****************************************************
  function norm_vec(a)

    type (vector), intent (in) :: a
    real(dp)                   :: norm_vec  

    norm_vec = norm2(a%x(:))

  end function norm_vec

  !****************************************
  ! returns the square of the given vector
  !****************************************
  function square(a)

    type(vector), intent (in) :: a
    real(dp)                  :: square

    square =    sum(a%x(:)**2)

  end function square

  !*****************************************************
  ! get the vector entries as array (can simply use a%x)
  !*****************************************************
  function get_array(a)

    type(vector), intent(in) :: a
    real(dp)                 :: get_array(NUM_SPAT_DIM)

    get_array=a%x

  end function get_array

  !*****************************************************
  ! get the entries of a VECTOR of VECTOR as an array
  !*****************************************************
  function get_array_1d(a)

    type(vector), dimension(:) :: a

    real(dp)                 :: get_array_1d(NUM_SPAT_DIM*size(a))
    integer                  :: j, is_j, ie_j, n

    n = size(a)

    do j = 1, n

       call split(j,is_j,ie_j) ! split j index storage

       get_array_1d(is_j:ie_j) = get_array(a(j))

    end do

  end function get_array_1d

  !*******************************************************************************************************
  ! constructor for a new matrix with entries supplied as an array (only NUM_SPAT_DIM compatible matrices)
  !*******************************************************************************************************
  function new_matrix_from_array(a)

    real(dp), intent(in) :: a(NUM_SPAT_DIM**2)
    type(matrix)         :: new_matrix_from_array

    new_matrix_from_array%ij = &
         &reshape(a, (/ NUM_SPAT_DIM, NUM_SPAT_DIM /))

  end function new_matrix_from_array

  !**********************************************************************************
  ! constructor for a new matrix with entries supplied as an array (any n by n matrix
  !**********************************************************************************
  function mat(arr,n)

    real(dp), intent(in)    :: arr(n*n)
    integer(sp), intent(in) :: n
    real(dp)                :: mat(n,n)

    mat   = reshape(arr, (/n,n/))

  end function mat

  ! get the matrix entries as array
  function get_matrix(A)
    type(matrix), intent(in) :: A
    real(dp)                 :: get_matrix(NUM_SPAT_DIM, NUM_SPAT_DIM)
    get_matrix = A%ij
  end function get_matrix

  ! unwraps a 2d matrix and stores as real numbers
  function get_matrix_2d(A) result(mat)

    type(matrix), dimension(:,:), intent(in) :: A

    real(dp)    :: mat(size(A,1)*NUM_SPAT_DIM, size(A,2)*NUM_SPAT_DIM) 
    integer(sp) :: i,j, m, n
    integer(sp) :: is_i, ie_i, is_j, ie_j

    n =size(A,1)
    m =size(A,2)

    do i = 1, n
       do j = 1, m

          call split(j,is_j,ie_j) ! split j index storage
          call split(i,is_i,ie_i) ! split i index storage

          mat(is_j:ie_j,is_i:ie_i) = get_matrix(A(j,i))

       end do
    end do

  end function get_matrix_2d

  ! unwraps a vector of vector and stores as an array
  function get_vector_elements(a,n)
    integer(sp) :: n
    type(vector), intent(in) :: a(n)
    real(dp)    :: get_vector_elements(n*NUM_SPAT_DIM) 
    integer(sp) :: i
    integer(sp) :: is_i, ie_i

    do i = 1, n
       call split(i,is_i,ie_i) ! split i index storage
       get_vector_elements(is_i:ie_i) = a(i)%x ! get_matrix(A(j,i))
       !        print*, is_j,ie_j,is_i,ie_i,i,j
    end do

  end function get_vector_elements

  subroutine split(i,is,ie)
    integer(sp) :: i
    integer(sp ):: is, ie
    is = (NUM_SPAT_DIM*(i-1) ) + 1
    ie = NUM_SPAT_DIM*( (i-1) + 1)
    !   print*, is,ie !,NUM_SPAT_DIM,i
  end subroutine split

  !------------------------
  ! transpose of a matrix
  !------------------------
  function trans(A)
    type(matrix),intent(in) :: A
    type(matrix)            :: trans
    trans = matrix(transpose(get_matrix(A)))
  end function trans

  !-------------------------------------------------------
  !returns the matrix addition of matrices of TYPE matrix
  !-------------------------------------------------------x
  function add_matrices(A, B)
    type(matrix), intent(in) :: A, B
    type(matrix)  :: add_matrices
    add_matrices%ij = A%ij +B%ij
  end function add_matrices

  !-------------------------------------------------------
  !returns the matrix addition of vectors of TYPE matrix
  !-------------------------------------------------------
  function add_vectors(a, b)
    type(vector), intent(in) :: a, b
    type(vector)  :: add_vectors
    add_vectors%x = a%x + b%x
  end function add_vectors

  !-------------------------------------------------------
  !returns the matrix subtraction of matrices of TYPE matrix
  !-------------------------------------------------------
  function sub_matrices(A, B)
    type(matrix), intent(in)  :: A, B
    type(matrix)  :: sub_matrices
    sub_matrices%ij = A%ij -B%ij
  end function sub_matrices

  !-------------------------------------------------------
  !returns the matrix addition of vectors of TYPE matrix
  !-------------------------------------------------------
  function sub_vectors(a, b)
    type(vector), intent(in) :: a, b
    type(vector)  :: sub_vectors
    sub_vectors%x = a%x - b%x
  end function sub_vectors

  !-------------------------------------------------------
  ! returns the negative of type matrix
  !-------------------------------------------------------
  function negate_matrix(A)
    type(matrix), intent(in)  :: A
    type(matrix)  :: negate_matrix
    negate_matrix%ij = -A%ij
  end function negate_matrix

  !-------------------------------------------------------
  ! returns the negative of type vector
  !-------------------------------------------------------
  function negate_vector(a)
    type(vector), intent(in)  :: a
    type(vector)  :: negate_vector
    negate_vector%x = -a%x
  end function negate_vector

  !-------------------------------------------------------
  ! convert from degree to radian
  !-------------------------------------------------------
  function deg2rad(deg)
    real(dp) :: deg2rad, deg
    deg2rad  =  deg*RAD_PER_DEG
  end function deg2rad

  !-------------------------------------------------------
  ! convert from radian to degree
  !-------------------------------------------------------
  function rad2deg(rad)
    real(dp) :: rad2deg, rad
    rad2deg  =  rad*DEG_PER_RAD
  end function rad2deg

  !************************************************
  ! returns an identity matrix of size NUM_SPAT_DIM
  !************************************************
  function idM()
    type(matrix) :: idM
    idM = matrix(eye(NUM_SPAT_DIM))
  end function idM

  !************************************************
  ! returns an zero matrix of size NUM_SPAT_DIM  
  !************************************************
  function zeroM()
    type(matrix) :: zeroM
    zeroM = matrix(zeros(NUM_SPAT_DIM))
  end function zeroM

  !************************************************
  ! returns an unit matrix of size NUM_SPAT_DIM
  !************************************************
  function unitM()
    type(matrix) :: unitM
    unitM =  matrix(ones(NUM_SPAT_DIM))
  end function unitM

  !************************************************
  ! generates a nxn identity matrix
  !************************************************
  function eye(n)
    integer(sp) :: n
    real(dp)    :: eye(n,n)
    integer     :: i, j
    eye = 0.0_dp
    do i = 1, n
       do j = 1, n
          if (j.eq.i) eye(j,i) = 1.0_dp
       end do
    end do
  end function eye

  !************************************************
  ! generates a nxn zero matrix
  !************************************************
  function zeros(n)
    integer(sp) :: n
    real(dp)    :: zeros(n,n)
    zeros        = 0.0_dp
  end function zeros

  !************************************************
  ! generates a nxn unit matrix
  !************************************************
  function ones(n)
    integer(sp) :: n
    real(dp)    :: ones(n,n)
    ones        = 1.0_dp
  end function ones

  

  !---------------------------------------------------------------------------
  ! function that reverses the order of the constituent elements in the vector
  !---------------------------------------------------------------------------
  function reverse_real(a) result(rev_a)
    real(dp), intent(in)     :: a(:)
    real(dp)                 :: rev_a(size(a))
    integer(sp)              :: i , n

    n= size(a)
    do i =  1, n
       rev_a(i) = a(n+1-i)
    end do

  end function reverse_real

  !---------------------------------------------------------------------------
  ! function that reverses the order of the constituent elements in the vector
  !---------------------------------------------------------------------------
  function reverse_int(a) result(rev_a)
    integer(sp), intent(in)     :: a(:)
    integer(sp)                 :: rev_a(size(a))
    integer(sp)                 :: i , n

    n= size(a)
    do i =  1, n
       rev_a(i) = a(n+1-i)
    end do

  end function reverse_int

!!$subroutine disp(t,y)
!!$  real(dp)    :: t, y(12)
!!$  integer(sp) :: i
!!$  write(*,'(13f13.2)') t, (y(i),i=1,12)
!!$end subroutine disp
!!$
!!$
!!$subroutine disp_mat(t, A,m,n)
!!$  implicit none
!!$  real(dp), intent(in) :: t
!!$  integer(sp)          :: i, j, m, n
!!$  real(dp)             :: A(m,n)
!!$  do i = 1, n
!!$     write(*,'(13f13.2)') t, (A(j,i),j=1,m)
!!$  end do
!!$end subroutine disp_mat

end module utils
