module constants
implicit none

!-----------------------------------------
! defining all the machine related values
!----------------------------------------

integer, parameter           :: dp = kind(0.d0)                   ! double precision
integer, parameter           :: sp = kind(0.0)                    ! single precision

!------------------------------------
! define all solver related constants
!------------------------------------
integer(sp), parameter       :: dim_q = 2                         ! number of variables in q vector q = [q1, q2, q3, ..., q_{dim_q}]
integer(sp), parameter       :: dim_x = 10                        ! number of variables in x_vector x = [x1, x2, x3, ..., x_{dim_x}]
real(dp), parameter          :: del_t = 0.1_dp                    ! time step used in integration

end module constants
