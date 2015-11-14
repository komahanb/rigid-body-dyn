!=====================================================================!
! Program that finds the first derivative of a function using complex 
! step method, finite-differencing and compares it with exact value.
!---------------------------------------------------------------------!
! Example:
!---------------------------------------------------------------------!
! Let f(z) = f(x+iy) = sin(z) then f'(x) = imag(f(x+ih))
!
! where "h" is simply a step size. The precision of the derivative is 
! EQUAL to the exponent of the step size. 
! 
! Since we know the  analytical  derivative of sin(x) to be cos(x), we
! can see the number of significant digits the derivatives would 
! agree to.
!=====================================================================!
program complex_step

  implicit none

  !-------------------------------------------------------------------!
  !  Define constants to manage precision [TUNABLE]
  !-------------------------------------------------------------------!

  integer, parameter      :: sp = kind(0.0)  ! single precision
  integer(sp), parameter  :: dp = kind(0.0d0)! double precision

  !-------------------------------------------------------------------!
  !  Problem parameters
  !-------------------------------------------------------------------!

  complex(dp)             :: z               ! define a complex number
  real(dp)                :: x, y, h         ! real, imag, step size

  ! define step size
  h = 1.0e-16_dp

  ! what is the real component
  x = 2.0_dp     

  ! MUST be zero for this method to work
  y = 0.0_dp      

  ! perturb the imaginary part with smallest possible val
  y = y + h       

  !make a complex number with real and perturbed imaginary part
  z = dcmplx(x,y) 

  write(*,*), "exact f'(x)", "                complex step", &
       &"             finite difference"

  print*, cos(x), imag(sin(z))/h, (sin(x+h)-sin(x))/h

end program complex_step
