program second_order
  use constants
  use solver_utils
  use dim_flexible_multi_body_dyn
  use system_components
  use utils
  implicit none

  real(dp), parameter       :: t0=0, t_final=1
  real(dp)                  :: current_time
  integer, parameter        :: nSteps=10 ! number of time steps
  integer                   :: i,j,k
  real(dp)                  :: q0(dim_q)
  real(dp)                  :: q(0:nSteps,dim_q)
  real(dp)                  :: q_dot(0:nSteps,dim_q)
  real(dp)                  :: q_double_dot(0:nSteps, dim_q)  ! state variable, its first and second derivatives
  real(dp), parameter       :: tol = 1.0d-8
  real(dp)                  :: residual
  real(dp)                  :: del_q(dim_q)
  
  type(body)                :: body_A
  type(body_fixed_frame)    :: frame_A
  type(matrix) :: test(1,1)

!  print*, size(test,1)
!  print*, get_matrix(test)
!  call split(4,i,j)
!  print*, get_matrix_2d(test,1,1)

  stop
!  call test_skew_sym
!  call test_cross_pdt
  call test_eye
  call test_ones
  call test_zeros
  call test_add_sub_mat
  
  stop
  

  frame_A%r_alpha     = (/ 1._dp, 1._dp, 1._dp /)
  frame_A%theta_alpha = (/ 30.0d0, 60.0d0, 0.0d0 /)/rad_to_deg
  frame_A%v_alpha     = (/ 1._dp, 1._dp, 1._dp /)
  frame_A%omega_alpha = (/ 1._dp, 1._dp, 1._dp /)
  frame_A%rot_mat     = rotate_to_body_frame(frame_A%theta_alpha)

  body_A%v_alpha   = matmul(frame_A%rot_mat, frame_A%v_alpha)
  print*, frame_A%v_alpha
  print*,  body_A%v_alpha
  print*,""

  
!!$
!!$  body_A%v_alpha   = matmul(frame_A%v_alpha, frame_A%rot_mat)
!!$  print*, frame_A%v_alpha
!!$  print*,  body_A%v_alpha
  stop
  !-------------------
  ! Initial conditions
  !-------------------

  q(0,:)     = (/ 1.0_dp, 1.0_dp /)           ! set initial condition for q(x=1.0) 
  q_dot(0,:) =  (/ 1.0_dp, 1.0_dp /)      ! set initial condition for q_dot
  
  ! (1)  get approximated q_double_dot based on the initial conditions
  !  q_double_dot(0,:) = get_approximated_q_double_dot(q(k-(m+d-1):k,:), m)
  q_double_dot(0,:) = (/ 1.0_dp, 1.0_dp/)
  
  del_q = 0.0_dp
  k = 0
  current_time = t0
 
  do while(current_time.eq.t_final)
     ! increase the time step and its counter
     current_time = current_time + del_t
     k = k + 1

     ! (2) form R Matrix (rhs) (-R)
     ! (3) form Jacobian (lhs) (partial_R/partial_q + (\alpha0/del_t) partial_R/partial_q_dot  + (\alpha0/del_t^2) partial_R/partial_q_double_dot)

     ! extrapolated q that is used as a good starting point for newton krylov
     q0 = get_extrapolated_q(q(k-1,:),q_dot(k-1,:),q_double_dot(k-1,:))

     !-----------------------------------------------
     ! (4) use a solver to solve the system for del_q
     !-----------------------------------------------
     !     del_q              = get_del_q()

     !-------------------------------------
     ! (5) update, q, q_dot_q_double_dot
     !-------------------------------------
     q(k,:)            = get_updated_q(q(k-1,:), del_q)
     q_dot(k,:)        = get_updated_q_dot(q_dot(k-1,:), del_q)
     q_double_dot(k,:) = get_updated_q_double_dot(q_double_dot(k-1,:), del_q)         
     
  end do
  
  print*, "time integration complete"
  
  
  ! adjoint pseudo code
  
  !=====================
  ! do unit testing here
  !=====================

  !call test_get_extrapolated_q
  !call test_get_bdf_coeffs
  call test_get_approximated_q_dot
  call test_get_approximated_q_double_dot

end program second_order

!-------------------------------------
! test second derivative approximation
!-------------------------------------
subroutine test_get_approximated_q_double_dot
  use constants
  use solver_utils
  implicit none

  real(dp)      :: q(0:4,dim_q)    ! I will simulate upto 5 available time-set data
  real(dp)      :: q_double_dot(dim_q)    ! stores the time-derivative
  integer(sp)   :: d, m, k

  ! first setup a q matrix with some simulated data
  q(0,1:dim_q) = (/ 0.1_dp, (1.1_dp)**3 /)
  q(1,1:dim_q) = (/ 0.2_dp, (1.2_dp)**3 /)
  q(2,1:dim_q) = (/ 0.3_dp, (1.3_dp)**3 /)
  q(3,1:dim_q) = (/ 0.4_dp, (1.4_dp)**3 /)
  q(4,1:dim_q) = (/ 0.5_dp, (1.5_dp)**3 /)
!!$

  q(0,1:dim_q) = (/ 0.1_dp, exp(0.1_dp) /)
  q(1,1:dim_q) = (/ 0.2_dp, exp(0.2_dp) /)
  q(2,1:dim_q) = (/ 0.3_dp, exp(0.3_dp) /)
  q(3,1:dim_q) = (/ 0.4_dp, exp(0.4_dp) /)
  q(4,1:dim_q) = (/ 0.5_dp, exp(0.5_dp) /)

  print*, "------------------------------------------------------"
  print*, "---------test_get_approximated_q_double_dot------------"
  print*, "------------------------------------------------------"  

  ! call the routine to and check the q_double_dot vector
  d=2; m=1; k=4
  print*, "d=2, m=1,k=4(t=0.4s)", get_approximated_q_double_dot(q(k-(m+d-1):k,:), m)

  d=2; m=2; k=4;
  print*, "d=2, m=2,k=4(t=0.4s)", get_approximated_q_double_dot(q(k-(m+d-1):k,:), m)

  d=2; m=3; k=4;
  print*, "d=2, m=3,k=4(t=0.4s)", get_approximated_q_double_dot(q(k-(m+d-1):k,:), m)


!  print*, "act  d=2,k=4(t=0.4s)", 0.0_dp, 6._dp *1.5_dp

  print*, "act  d=2,k=4(t=0.4s)", 0.0_dp, exp(0.5_dp)

  !  d=1; m=3
  !  print*, "d=1, m=3", get_approximated_q_double_dot(q(0:d+m-1,:), 5 , m)
!!$  d=1; m=4
!!$  print*, "d=1, m=4", get_approximated_q_double_dot(q(0:d+m-1,:), 5 , m)
  print*, "------------------------------------------------------"
  print*, ""
end subroutine test_get_approximated_q_double_dot

!-------------------------------------
! test first derivative approximation
!-------------------------------------
subroutine test_get_approximated_q_dot
  use constants
  use solver_utils
  implicit none

  real(dp)      :: q(0:4,dim_q)    ! I will simulate upto 5 available time-set data
  real(dp)      :: q_dot(dim_q)    ! stores the time-derivative
  integer(sp)   :: d, m, k

  ! first setup a q matrix with some simulated data
  q(0,1:dim_q) = (/ 0.1_dp, (1.1_dp)**3 /)
  q(1,1:dim_q) = (/ 0.2_dp, (1.2_dp)**3 /)
  q(2,1:dim_q) = (/ 0.3_dp, (1.3_dp)**3 /)
  q(3,1:dim_q) = (/ 0.4_dp, (1.4_dp)**3 /)
  q(4,1:dim_q) = (/ 0.5_dp, (1.5_dp)**3 /)


  q(0,1:dim_q) = (/ 0.1_dp, exp(0.1_dp) /)
  q(1,1:dim_q) = (/ 0.2_dp, exp(0.2_dp) /)
  q(2,1:dim_q) = (/ 0.3_dp, exp(0.3_dp) /)
  q(3,1:dim_q) = (/ 0.4_dp, exp(0.4_dp) /)
  q(4,1:dim_q) = (/ 0.5_dp, exp(0.5_dp) /)

  print*, "------------------------------------------------------"
  print*, "---------test_get_approximated_q_dot------------------"
  print*, "------------------------------------------------------"  

  ! call the routine to and check the q_dot vector
  d=1; m=1; k=4
  print*, "d=1, m=1,k=4(t=0.4s)", get_approximated_q_dot(q(k-(m+d-1):k,:), m)

  d=1; m=2; k=4;
  print*, "d=1, m=2,k=4(t=0.4s)", get_approximated_q_dot(q(k-(m+d-1):k,:), m)

  d=1; m=3; k=4;
  print*, "d=1, m=3,k=4(t=0.4s)", get_approximated_q_dot(q(k-(m+d-1):k,:), m)

  d=1; m=4; k=4;
  print*, "d=1, m=3,k=4(t=0.4s)", get_approximated_q_dot(q(k-(m+d-1):k,:), m)

!  print*, "act  d=1,k=4(t=0.4s)", 1.0_dp, 3._dp *1.5_dp*1.5_dp

  print*, "act  d=1,k=4(t=0.4s)", 1.0_dp, exp(0.5_dp)

  !  d=1; m=3
  !  print*, "d=1, m=3", get_approximated_q_dot(q(0:d+m-1,:), 5 , m)
!!$  d=1; m=4
!!$  print*, "d=1, m=4", get_approximated_q_dot(q(0:d+m-1,:), 5 , m)
  print*, "------------------------------------------------------"
  print*, ""
end subroutine test_get_approximated_q_dot

!----------------------------------------------------------------------
! routine that will check the coefficients need for 1 ,2 nd der approx
!----------------------------------------------------------------------
subroutine test_get_bdf_coeffs
  use constants
  use solver_utils
  implicit none

  print*, "------------------------------------------------------"
  print*, "---------test_get_bdf_coeffs--------------------------"
  print*, "------------------------------------------------------"
  print*, "case 1: derivative d = 1 and order m = 1"
  print*, get_bdf_coeffs(1, 1)
  print*, "case 2: derivative d = 1 and order m = 2"
  print*, get_bdf_coeffs(1, 2)
  print*, "case 3: derivative d = 1 and order m = 3"
  print*, get_bdf_coeffs(1, 3)
  print*, "------------------------------------------------------"
  print*, "case 4: derivative d = 2 and order m = 1"
  print*, get_bdf_coeffs(2, 1)
  print*, "case 5: derivative d = 2 and order m = 2"
  print*, get_bdf_coeffs(2, 2)
  print*, "case 6: derivative d = 2 and order m = 3"
  print*, get_bdf_coeffs(2, 3)
  print*, "case 7: derivative d = 2 and order m = 4"
  print*, get_bdf_coeffs(2, 4)
  print*, "--------completed--------------------------------------"
  print*, ""

end subroutine test_get_bdf_coeffs

!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
subroutine test
  use constants
  use solver_utils
  implicit none

end subroutine test

!----------------------------------------------------------------------
! routine to test the extrapolation of the first and second derivatives
!----------------------------------------------------------------------
subroutine test_get_extrapolated_q
  use constants
  use solver_utils
  implicit none

  real(dp):: q(2)
  real(dp) ::q0(2),q0_dot(2), q0_double_dot(2)
  real(dp):: x

  print*, "------------------------------------------------------"
  print*, "---------test_get_extrapolated_q----------------------"
  print*, "------------------------------------------------------"

  ! set constants
  x=1._dp

  ! give starting values for extrapolation
  q0=(/x**3, sin(x)/)
  q0_dot = (/3._dp*x**2, cos(x)/)
  q0_double_dot = (/6._dp*x, -sin(x)/)

  ! call the extrapolation routine
  print*, "old q=", q0
  q(:) = get_extrapolated_q( q0, q0_dot, q0_double_dot)
  print*, "new q=", q

  ! call the extrapolation routine
  print*, "old q=", q0
  q(:) = get_extrapolated_q( q0, q0_dot)
  print*, "new q=", q


  x=1._dp+del_t
  q0=(/ x**3, sin(x) /)

  print*, "ext q=", q0
  print*, "--------------complete--------------------------------"
  print*, ""

  return
end subroutine test_get_extrapolated_q

subroutine test_skew_sym
  use constants
  use utils
  implicit none

  real(dp)                  :: a(3)
  real(dp)                  :: a_skew(3,3)

  a = (/ 1, 2, 3/)
  a_skew = get_matrix(skew(vector(a)))

  print*, "a=",a
  print*, "anti-symmetric matrix:"
  print*, a_skew

end subroutine test_skew_sym

subroutine test_cross_pdt
  use constants
  use utils
  implicit none

  type(vector) ::a1, b1
  type(matrix) :: A
  real(dp) :: con =2

  A%ij(1,1) = 1
  A%ij(1,2) = 1
  A%ij(1,3) = 1
  A%ij(2,1) = 1
  A%ij(2,2) = 1
  A%ij(2,3) = 1
  A%ij(3,1) = 1
  A%ij(3,2) = 1
  A%ij(3,3) = 1

!  A%ij = con*A%ij

 ! print*, con*A%ij

!  stop

     a1%x = (/  -4,1,0 /)
     b1%x = (/ 4,1,0 /)
     
     !  a1=vector((/1,2,3/))
     !  b1=vector((/2,3,4/))
     !  real(dp) :: cross_pdt(3,3)

  print*, "a          =", a1
  print*, "b          =", b1
  print*, "cross pdt  =", cross(a1,b1)
  print*, "cross pdt  =", a1*b1
  print*, "Dot product=", dot(a1,b1)

end subroutine test_cross_pdt

subroutine test_eye
use constants
use matrix_utils
print*,"Testing identitiy matrix"

print*,"1x1:",eye(1)
print*,""

print*,"2x2:",eye(2)
print*,""

print*,"3x3:",matrix(eye(3))
print*,""


print*,"4x4:",eye(4)
print*,""

end subroutine test_eye


! test unit matrix
subroutine test_ones
use matrix_utils
print*,"Testing UNIT matrix"

print*,"1x1:",ones(1)
print*,""

print*,"2x2:",ones(2)
print*,""

print*,"3x3:",matrix(ones(3))
print*,""


print*,"4x4:",ones(4)
print*,""

end subroutine test_ones


subroutine test_zeros
use matrix_utils
print*,"Testing ZERO matrix"

print*,"1x1:",zeros(1)
print*,""

print*,"2x2:",zeros(2)
print*,""

print*,"3x3:",matrix(zeros(3))
print*,""


print*,"4x4:",zeros(4)
print*,""

end subroutine test_zeros



subroutine test_add_sub_mat
use constants
use utils
  real(dp) :: a1(3,3), b1(3,3)

  a1 = 0.0_dp
  b1 = -1.0_dp

  print* ,  a1-b1
  print* ,  matrix(a1)-matrix(b1)
end subroutine test_add_sub_mat
