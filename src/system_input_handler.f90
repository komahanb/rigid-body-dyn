!=====================================================================!
! Module to handle the inputs of the rigid-body-dynamics program
!---------------------------------------------------------------------!
! (a) setting up bodies, its state, its inetial props
! (b) setting up the number of bodies,
! (c) setting up the joints and the bodies they are connected to etc
!---------------------------------------------------------------------!
! Note:
!!--------------------------------------------------------------------!
! All the preprocessing must be done here and everything must be set 
! into the SYSYEM object (dynsys) referenced from global_variables.f90. 
! This dynsys object is used for further jacobian and residual assembly
! in jacobian.f90 and residual.f90
!=====================================================================!

module system_input_handler

  use global_constants, only: dp, sp, ZERO, TOT_NDOF, NUM_SPAT_DIM
  use global_variables, only: dynsys
  use rigid_body_class, only: rigid_body, print_rigid_body
  use joint_class, only: joint
  use pendulum_class, only : pendulum
  use dispmodule, only: disp

  implicit none
  private
  public :: read_system_input

contains 

  !*******************************************************************!
  ! Routine that reads the user input and sets the necessary variables
  !*******************************************************************!
  subroutine read_system_input()

    ! reference point on the body measured in body frame
    real(dp), dimension(NUM_SPAT_DIM)  :: c, fr, gr    

    real(dp), dimension(NUM_SPAT_DIM, NUM_SPAT_DIM)  :: J

    ! mass of the body
    real(dp)                           :: mass 

    ! the state
    real(dp)                           :: q(TOT_NDOF)= 0.0_dp

    ! time derivative of state
    real(dp)                           :: qdot(TOT_NDOF) = 0.0_dp

    type(rigid_body)                   :: body


    !-----------------------------------------------------------------!
    !----------------CREATE BODIES AND JOINTS-------------------------!
    !-----------------------------------------------------------------!

    !-----------------------------------------------------------------!
    ! read input from file or from an external program
    !-----------------------------------------------------------------!
    
    ! currently setting the values here directly

    call random_seed()
    call random_number(q)
    call random_number(qdot)
    !qdot = qdot**2

    ! set the body properties

    mass = 2.0_dp

    c= (/ 0.2, 0.3, 0.4 /)

    J = reshape((/ &
         &0.2, 0.3, 0.4, &
         &0.2, 0.3, 0.4, &
         &0.2, 0.3, 0.4 /), &
         &(/ 3,3 /))

    fr = (/ ZERO, ZERO, ZERO/)

    gr = (/ ZERO, ZERO, ZERO/)

    call disp(" >> Creating a body...")

    body = rigid_body(mass, c, J, fr, gr, q, qdot)

    call print_rigid_body(body)

    !-----------------------------------------------------------------!
    !----------------------- CREATE SYSTEM ---------------------------!
    !-----------------------------------------------------------------!

    call disp(" >> Creating a system...")
    
    ! Create a system by allocating the dynamic system as pendulum
    
    allocate(dynsys, source = pendulum(nbody = 1, njoint = 0))
    
    ! The above step just allocated space. Now we add the bodies and 
    ! joints into the system

    call dynsys % add_body(bnum = 1, bdy = body)

    call disp(" >> System is created...")

  end subroutine read_system_input

end module system_input_handler
