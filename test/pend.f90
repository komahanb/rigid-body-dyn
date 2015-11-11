!===========================================================
! program to run dynamics simulation of a physical pendulum
!===========================================================
program pendulum

  ! import the necessary modules
  !  use iso_c_binding ! to use tripan
  use global_constants
  use global_variables
  use types
  use utils
  use solver_utils
  use bodies, only: create_body, set_state
  use residual, only: get_residual
  use jacobian, only: get_jacobian
  use finite_diff
  implicit none

  ! position of the point of interest on the body measured in body frame
  real(dp), dimension(NUM_SPAT_DIM)  :: re    

  ! mass of the body
  real(dp)                           :: mass 

  ! the state
  real(dp)                           :: q(TOT_NDOF)= 0.0_dp
  real(dp)                           :: qtmp(TOT_NDOF)= 0.0_dp

  ! time derivative of state
  real(dp)                           :: q_dot(TOT_NDOF) = 0.0_dp
  real(dp)                           :: q_dot_tmp(TOT_NDOF) = 0.0_dp

  ! residual
  real(dp)                           :: res(TOT_NDOF)

  real(dp)                           :: res2(TOT_NDOF)

  ! jacobian
  real(dp)                           :: jac(TOT_NDOF, TOT_NDOF)
  real(dp)                           :: jac2(TOT_NDOF, TOT_NDOF)

  ! update from newton iterations
  real(dp)                           :: dq(TOT_NDOF) = 0.0_dp

  ! body -- single pendulum
  type(body)                         :: body1

  ! time, and norms to stop integration
  real(dp)                           :: time,l2_norm, linf_norm, dh

  ! loop variable
  integer(sp)                        :: k


  call disp("========================================================")
  call disp("                 Rigid body dynamics                    ")
  call disp("========================================================")

  !-------------------------------------------------------------------!
  ! Init routine that has sanity check on the input settings
  !-------------------------------------------------------------------!
  ! Note: 
  !
  ! Has nothing as of now, we can use to:
  ! Load input files, external calls etc
  !-------------------------------------------------------------------!
  call init() 

  !-------------------------------------------------------------------!
  ! (1) set the initial states and attributes of the body
  !-------------------------------------------------------------------!
  ! configure random number generator
  call random_seed()

  ! guess values
  call random_number(q)
  call random_number(q_dot)

  ! define the attributes of the body
  mass      = 2.0_dp
  ! point of interest on the body
  CALL RANDOM_NUMBER(re)
  !  re        = (/ 0.0_dp, 0.0_dp, 0.0_dp /)   

  !-------------------------------------------------------------------!
  ! (2) create a pendulum body using the state and attr above
  !-------------------------------------------------------------------!

  !  call disp(" >> Creating a body...")
  body1 = create_body(mass, vector(re), q, q_dot)
  !call disp(" >> Body created successfully...")

  !-------------------------------------------------------------------!
  ! Settings for time-integration
  !-------------------------------------------------------------------!

!!$! The following defaults are already set in global_varaibles.f90.
!!$! The user is free to change here too.
!!$  dT         = 0.01_dp
!!$  start_time = 0.0_dp
!!$  end_time   = 1.0_dp

  aa   = 1.0_dp/dT    ! used in jacobian assembly and state update
  bb   = 1.0_dp/dt**2 ! used in jacobian assembly and state update

  time = start_time

  march: do while ( time .le. end_time)  ! loop for time-marching

     time = time + dT ! update the time

     ! Newton iterations
     newton: do k = 1, MAX_NEWTON_ITER

        ! (3) Assemble the residual and jacobian using the body

        ! call this method to update the state for every other iter
        ! we don't have to create a whole new body again
        if (k .gt. 1) then
           call set_state(q, q_dot, body1)
        end if

        ! residual assembly
        !call disp(" >> Assembling the residuals...")
        res  = get_residual(body1)
        !call disp(" >> Residual assembly complete...")
        !call disp("   R   =   ", res)

        ! jacobian assembly
        !call disp(" >> Assembling the Jacobian...")
        jac  = get_jacobian(body1)
        !call disp(" >> Jacobian assembly complete...")
        !call disp("   JAC =   ",jac)

        ! ******************************************************
        ! (3) Solve the linear system and compute update delta_q
        ! ******************************************************

        !call disp(" >> Calling the linear solver...") 

        !dq = linear_solve(jac,-res,'GMRES')

        !call disp(" >> Update to residual is computed...") 

        l2_norm   = norm2(dq); linf_norm =  maxval(abs(dq));

        ! call disp(" >> L2-Norm of the update    :", l2_norm) 
        ! call disp(" >> Infty-Norm of the update :", linf_norm)

        !-------------------------------------------------------------!
        ! (4) Update the state variables and then update the body
        !-------------------------------------------------------------!

        q            = get_updated_q(q, dq)
        q_dot        = get_updated_q_dot(q_dot, dq)

        ! q_double_dot = get_updated_q_double_dot(q_double_dot, dq)

        ! call disp("   q      =   ", q)
        ! call disp("   qdot   =   ", q_dot)

        if (l2_norm .le. REL_TOL .OR. linf_norm.le. ABS_TOL) &
             & exit newton

        if (k .eq. MAX_NEWTON_ITER) then
           call disp(" >> Newton solution failed in ", k , "iterations")
           stop
        end if

     end do newton

     ! print the summary

  end do march

  call disp(" >> End of program...") 
  call disp("========================================================")

end program pendulum

!*********************************************************************!
! Callbaclk function for approximation of jacobian using finite-diff
!*********************************************************************!
subroutine residual2(q, qdot, f)

  use global_constants, only: sp, dp, TOT_NDOF, NUM_SPAT_DIM
  use types
  use bodies
  use residual

  real(dp), intent(in), dimension(TOT_NDOF)  :: q, qdot
  real(dp), intent(out), dimension(TOT_NDOF) :: f
  real(dp), dimension(NUM_SPAT_DIM)  :: re    
  real(dp)   :: mass
  type(body) :: body1

  mass  = 2.0_dp;   re   = (/ 0.1_dp, 0.2_dp, 0.3_dp /);
  body1 = create_body(mass, vector(re), q, qdot);

  f = get_residual(body1);  

!!$  f = (/ sin(q(1))*sin(q(3)), sin(q(2)), sin(q(3)), sin(q(4)),&
!!$       & sin(q(5)), sin(q(6)), &
!!$       &sin(q(7)), sin(q(8)), sin(q(9)), &
!!$       &sin(q(10)), sin(q(11)), sin(q(12)*sin(q(1))) /)
!!$

end subroutine residual2

!*********************************************************************!
! Callbaclk function for approximation of jacobian using finite-diff
!*********************************************************************!
subroutine residual(q, f)

  use global_constants, only: sp, dp, TOT_NDOF

  real(dp), intent(in), dimension(TOT_NDOF)  :: q
  real(dp), intent(out), dimension(TOT_NDOF) :: f

  !  f = get_residual();  

  f = (/ sin(q(1))*sin(q(3)), sin(q(2)), sin(q(3)), sin(q(4)),&
       & sin(q(5)), sin(q(6)), &
       &sin(q(7)), sin(q(8)), sin(q(9)), &
       &sin(q(10)), sin(q(11)), sin(q(12)*sin(q(1))) /)

end subroutine residual

!!$
!!$
!!$
!!$ 
!!$  ! residual assembly
!!$  !  call disp(" >> Assembling the residuals...")
!!$  res  = get_residual(body1)
!!$  !  call disp(" >> Residual assembly complete...")
!!$  !call disp("   R   =   ", res)
!!$
!!$  ! jacobian assembly
!!$  !call disp(" >> Assembling the Jacobian...")
!!$  jac  = get_jacobian(body1)
!!$  !call disp(" >> Jacobian assembly complete...")
!!$  !call disp("   JAC =   ",jac)
!!$
!!$
!!$  ! implement FD to check if the jacobian is the same
!!$  qtmp = q
!!$  call set_state(q, q_dot, body1)
!!$  dh = 1.0d-8
!!$
!!$  res = get_residual(body1)
!!$     
!!$  var: do  k = 1 , TOT_NDOF
!!$
!!$     ! perturb k-th variable
!!$     qtmp( k ) = qtmp (k)  + dh
!!$     ! updating the body state
!!$     call set_state(qtmp, q_dot, body1)
!!$     ! find the solution
!!$     res2 = get_residual(body1)
!!$     !    call disp("   R   =   ", res)
!!$     
!!$     dq =  (res2-res)/dh
!!$     
!!$     jac2 (:, k) = dq
!!$     
!!$     qtmp = q
!!$
!!$  end do var
!!$  call disp("   dJAC =   ",jac2)
!!$ 
!!$  call disp("   Error =   ",jac2-Jac)
!!$
!!$
!!$
!!$  stop"stopped"
!!$  ! define the initial states
!!$  q(1:3)    = (/ 1.0_dp, 2.0_dp, 3.0_dp /)
!!$  q(4:6)    = (/ deg2rad(20.0_dp), deg2rad(0.0_dp), deg2rad(20.0_dp) /)
!!$  q(7:9)    = (/ 1.0_dp,  2.0_dp, 1.0_dp /)
!!$  q(10:12)  = (/ 2.0_dp, 4.0_dp, 5.0_dp /)
!!$
!!$  ! define the intial time derivative of state
!!$  q_dot(1:3)   = (/ 2.0_dp, 2.0_dp, 2.0_dp /)
!!$  q_dot(4:6)   = (/ 0.25_dp, 0.0_dp, 0.0_dp /)
!!$  q_dot(7:9)   = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
!!$  q_dot(10:12) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)


subroutine matset_triad ( n, mode, nelt, ia, ja, a )

  !*****************************************************************************80
  !
  !! MATSET_TRIAD sets the data structure for a sparse matrix.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 July 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer N, the number of elements in the vectors.
  !
  !    Input, integer MODE.
  !    0, setup mode.  Count the number of nonzero elements.
  !    1, assignment mode.  NELT is set.  Set the matrix entries.
  !
  !    Input/output, integer NELT, the number of nonzero entries.
  !    On call with MODE = 0, this is an output quantity.
  !    On call with MODE = 1, this is an input quantity.
  !
  !    Output, integer IA(NELT), JA(NELT), real ( kind = 8 ) A(NELT), the data
  !    structure storing the sparse matrix.  These values are only
  !    set on a call with MODE = 1.
  !
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nelt

  real ( kind = 8 ) a(nelt)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(nelt)
  integer ( kind = 4 ) ja(nelt)
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) k

  k = 0
  do i = 1, n

     if ( 1 < i ) then
        k = k + 1
        if ( mode == 1 ) then
           ia(k) = i
           ja(k) = i-1
           a(k) = -1.0D+00
        end if
     end if

     k = k + 1
     if ( mode == 1 ) then
        ia(k) = i
        ja(k) = i
        a(k) = 2.0D+00
     end if

     if ( i < n ) then
        k = k + 1
        if ( mode == 1 ) then
           ia(k) = i
           ja(k) = i+1
           a(k) = -1.0D+00
        end if
     end if

  end do

  if ( mode == 0 ) then
     nelt = k
  end if

  return
end subroutine matset_triad
subroutine msolve_identity ( n, r, z, nelt, ia, ja, a, isym, rwork, iwork )

  !*****************************************************************************80
  !
  !! MSOLVE_IDENTITY applies the identity matrix preconditioner.
  !
  !  Discussion:
  !
  !    Most SLAP solver routines require a preconditioner routine
  !    that can solve M * Z = R.  If no preconditioning is required,
  !    then you can simply behave as though the preconditioning matrix
  !    M was the identity matrix.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    21 July 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of elements in the vectors.
  !
  !    Input, real ( kind = 8 ) R(N), the right hand side.
  !
  !    Output, real ( kind = 8 ) Z(N), the solution of M * Z = R.
  !
  !    Input, integer ( kind = 4 ) NELT, the number of nonzero entries in A.
  !
  !    Input, integer ( kind = 4 ) IA(NELT), JA(NELT), real ( kind = 8 ) A(NELT), 
  !    the data structure storing the sparse matrix.
  !
  !    Input, integer ( kind = 4 ) ISYM, is 0 if all nonzero entries of the matrix
  !    are stored, and 1 if only the diagonal and upper or lower triangle
  !    are stored.
  !
  !    Input, real ( kind = 8 ) RWORK(*), a real array that
  !    can be used to pass information to the preconditioner.
  !
  !    Input, integer ( kind = 4 ) IWORK(*), an integer array that
  !    can be used to pass information to the preconditioner.
  !
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nelt

  real ( kind = 8 ) a(nelt)
  integer ( kind = 4 ) ia(nelt)
  integer ( kind = 4 ) isym
  integer ( kind = 4 ) iwork(*)
  integer ( kind = 4 ) ja(nelt)
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) rwork(*)
  real ( kind = 8 ) z(n)

  z(1:n) = r(1:n)

  return
end subroutine msolve_identity
