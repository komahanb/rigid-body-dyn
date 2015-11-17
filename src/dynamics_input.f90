!=====================================================================!
! Module to handle the inputs of the rigid-body-dynamics program
! i.e. setting up bodies, its state, its inetial props
! setting up the number of bodies,
! setting up the joints and the bodies they are connected to etc...
!=====================================================================!

module dynamics_input

  use global_constants
  use global_variables
  use body_class
  use rigid_body_class
  use joint_class
  use types
  
  implicit none
  private
  public :: read_input 

  ! residual
contains 

  !*******************************************************************!
  ! Routine that reads the user input and sets the necessary variables
  !*******************************************************************!
  subroutine read_input()

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
    ! read input from file or from an external program


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
    
    call set_body_num(body, 2)
    call set_body_type(body, "CROD1")

    print*, get_body_num(body)
    print*, get_body_type(body)


    call print_rigid_body(body)

    

  end subroutine read_input

end module dynamics_input
