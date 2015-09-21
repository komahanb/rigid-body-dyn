module test
  use utils

  ! each body shall contain the following
  ! body fixed frame located at O position(R, Theta) wrt inertial
  ! the velocity v0 (body)
  ! the angular omega (body)
  ! A rotation matrix (inertia-body)
  ! A transformation matrix (inertial-body)

  ! mass moment of inertia
  ! stiffness matrix
  ! EQ1, EQ2, EQ3 = m{v0_dot} - {c} x {omega_dot}  + {p}.{qs_double_dot} + {omega} x m {v0} - {omega} x {c} x {omega} + {omega} x {p}.{qs_dot}
  ! EQ4, EQ5, EQ6 = {c}.{v0_dot}  + [J]{omega_dot} + [h]{qs_double_dot} + [c] x {omega} x {v0} + {omega} x [J] {omega} + {}x[p].{qs_dot} + {omega} x [h]{qs_dot} =gr

  type body
 
     type(vector) :: r
     type(vector) :: theta

     type(vector) :: v, v_dot
     type(vector) :: omega, omega_dot

     real(dp)     :: m, C(num_spat_dim, num_spat_dim)

     type(vector) :: qs, qs_dot, qs_double_dot

     real(dp)     :: D(num_spat_dim, num_spat_dim)
     real(dp)     :: S(num_spat_dim, num_spat_dim)

     real(dp)     :: J(num_spat_dim, num_spat_dim)
     real(dp)     :: K(num_spat_dim, num_spat_dim)
     real(dp)     :: P(num_spat_dim, num_spat_dim)
      
  end type body
  
end module test
