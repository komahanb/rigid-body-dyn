module system_components
  use constants

  type body_fixed_frame 
     ! location
     real(dp), dimension(num_spat_dim)  :: r_alpha      ! position of the frame wrt inertial
     real(dp), dimension(num_spat_dim)  :: theta_alpha  ! rotation of the frame wrt inertial
     ! time derivatives
     real(dp), dimension(num_spat_dim)  :: v_alpha      ! d/dt (r_alpha) wrt inertial
     real(dp), dimension(num_spat_dim)  :: omega_alpha  ! d/dt (theta_alpha) wrt inertial
     real(dp), dimension(num_spat_dim, num_spat_dim) :: rot_mat
  end type body_fixed_frame

  ! define a body
  type body
     type(body_fixed_frame)  :: frame        ! each body has a body-fixed ref frame
     real(dp), dimension(num_spat_dim)                :: v_alpha      ! velocity of body wrt body-fixed frame
     real(dp), dimension(num_spat_dim)                :: omega_alpha  ! angular velocity of body wrt body-fixed frame
     real(dp), dimension(num_spat_dim)                :: f_alpha_beta ! reaction force on alpha from beta wrt body-fixed frame
     real(dp), dimension(num_spat_dim)                :: g_alpha_beta ! reaction torque on alpha from beta wrt body-fixed frame
     real(dp), dimension(num_spat_dim)                :: r_alpha_i    ! position vector of the connection joint from the body frame
     real(dp)                :: m_alpha      ! mass
     real(dp)                :: c_alpha
     real(dp)                :: J_alpha
     real(dp)                :: q_r_alpha
     real(dp)                :: q_s_alpha
  end type body

contains

  ! rotates a vector from inertial to body frame
  function rotate_to_body_frame(theta) result(c_mat)
    
    real(dp), intent(in)    :: theta(num_spat_dim)
    real(dp)                :: c_mat(num_spat_dim, num_spat_dim)

    c_mat(1,1) =  cos(theta(2))*cos(theta(3)) + sin(theta(1))*sin(theta(2))*sin(theta(3))
    c_mat(1,2) = -cos(theta(2))*sin(theta(3)) + sin(theta(1))*sin(theta(2))*cos(theta(3))
    c_mat(1,3) =  cos(theta(1))*sin(theta(2))

    c_mat(2,1) =  cos(theta(1))*sin(theta(3))
    c_mat(2,2) =  cos(theta(1))*cos(theta(3))
    c_mat(2,3) = -sin(theta(1))

    c_mat(3,1) = -sin(theta(2))*cos(theta(3)) + sin(theta(1))*cos(theta(2))*sin(theta(3))
    c_mat(3,2) =  sin(theta(2))*sin(theta(3)) + sin(theta(1))*cos(theta(2))*cos(theta(3))
    c_mat(3,3) =  cos(theta(1))*cos(theta(2))

  end function rotate_to_body_frame

end module system_components
