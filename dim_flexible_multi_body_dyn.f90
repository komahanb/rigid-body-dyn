module dim_flexible_multi_body_dyn
  implicit none

  integer, parameter:: print_level=3

contains

  ! returns alpha_i for the order
  function get_alpha(order) result(alpha)
    implicit none

    real(kind=8)::alpha
    integer::order

    alpha = dble(order)

  end function get_alpha

  ! returns beta_i for the order
  function get_beta(order) result(beta)
    implicit none

    real(kind=8)::beta
    integer::order

    beta = dble(order)

  end function get_beta


end module dim_flexible_multi_body_dyn
