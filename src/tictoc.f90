!=====================================================================!
! Module represents a timer that can be started and stopped to measure
! how long different routines run. The intrinsic routine system_clock
! is used to measure time rather than cpu_time.
!=====================================================================!
module tictoc

  use global_constants, only: ZERO, dp

  implicit none

  type timer
     logical  :: running      = .false. ! is timer running?
     integer  :: start_counts = 0       ! counts when started
     real(dp) :: elapsed      = 0.0_dp  ! total time elapsed in seconds
!!$   contains
!!$     procedure :: start     => timer_start
!!$     procedure :: get_value => timer_get_value
!!$     procedure :: stop      => timer_stop
!!$     procedure :: reset     => timer_reset
  end type timer

contains

  !*******************************************************************!
  ! TIMER_START
  !*******************************************************************!

  subroutine timer_start(self)

    type(timer), intent(inout) :: self

    ! Turn timer on and measure starting time

    self % running = .true.

    call system_clock(self % start_counts)

  end subroutine timer_start

  !*******************************************************************!
  ! TIMER_GET_VALUE
  !*******************************************************************!

  function timer_get_value(self) result(elapsed)

    type(timer), intent(in) :: self   ! the timer
    real(8)                 :: elapsed ! total elapsed time

    integer :: end_counts   ! current number of counts
    integer :: count_rate   ! system-dependent counting rate
    real    :: elapsed_time ! elapsed time since last start

    if (self % running) then

       call system_clock(end_counts, count_rate)

       elapsed_time = real(end_counts &
            & - self % start_counts)/real(count_rate)

       elapsed = self % elapsed + elapsed_time

    else

       elapsed = self % elapsed

    end if

  end function timer_get_value


  !*******************************************************************!
  ! TIMER_STOP
  !*******************************************************************!

  subroutine timer_stop(self)

    type(timer), intent(inout) :: self

    ! Check to make sure timer was running
    if (.not. self % running) return

    ! Stop timer and add time

    self % elapsed = timer_get_value(self)

    self % running = .false.

  end subroutine timer_stop

  !*******************************************************************!
  ! TIMER_RESET
  !*******************************************************************!

  subroutine timer_reset(self)

    type(timer), intent(inout) :: self

    self % running      = .false.
    self % start_counts = 0
    self % elapsed      = ZERO

  end subroutine timer_reset

end module tictoc
