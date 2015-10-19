module constants
implicit none

!-----------------------------------------
! defining all the machine related values
!----------------------------------------

integer, parameter           :: dp = kind(0.d0)                   ! double precision
integer, parameter           :: sp = kind(0.0)                    ! single precision
integer, parameter           :: units_len = 10                    ! length of the units

integer, parameter           :: ref_frame_len = 8

integer, parameter           :: num_spat_dim  = 3 

real(dp), parameter          :: pi = 22._dp/7._dp
real(dp), parameter          :: rad_to_deg = 180._dp/pi

real(dp), parameter          :: dzero=0._dp
real(dp), parameter          :: done=1._dp




!------------------------------------
! define all solver related constants
!------------------------------------
integer(sp), parameter       :: dim_q = 2                         ! number of variables in q vector q = [q1, q2, q3, ..., q_{dim_q}]
integer(sp), parameter       :: dim_x = 10                        ! number of variables in x_vector x = [x1, x2, x3, ..., x_{dim_x}]
real(dp), parameter          :: del_t = 0.1_dp                    ! time step used in integration

type vector
   real(dp)                  :: x(num_spat_dim) = 0.0_dp
end type vector

type matrix
   real(dp)                  :: ij(num_spat_dim, num_spat_dim) = 0.0_dp
end type matrix

end module constants
