module dynamics_input

  implicit none
  private
  public  :: read_input 

  ! reference point on the body measured in body frame
  real(dp), dimension(NUM_SPAT_DIM)  :: re    

  ! mass of the body
  real(dp)                           :: mass 

  ! the state
  real(dp)                           :: q(TOT_NDOF)= 0.0_dp

  ! time derivative of state
  real(dp)                           :: q_dot(TOT_NDOF) = 0.0_dp

  ! residual
contains 

  !*******************************************************************!
  ! Routine that reads the user input and sets the necessary variables
  !*******************************************************************!
  subroutine read_input()

    ! read input from file or from an external program


    ! currently setting the values here directly

    call random_seed()
    call random_number(q)
    call random_number(qdot)
    !qdot = qdot**2

    ! mass of the body
    mass = 2.0_dp

    ! used to calculate the inertial properties J and C
    re   = (/ 0.1_dp, 0.2_dp, 0.3_dp /) 

    call disp(" >> Creating a body...")

    body1 = create_body(mass, vector(re), q, q_dot)

    call print_body(body1)

  end subroutine read_input

end module dynamics_input
