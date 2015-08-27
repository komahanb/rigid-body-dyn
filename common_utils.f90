module utils
  implicit none
!  private
!  public newunit

contains

!------------------------------------------------------------
! a function that returns a new unit number for file handling
!------------------------------------------------------------
  integer function newunit(unit) result(n)
    ! returns lowest i/o unit number not in use
    integer, intent(out), optional :: unit
    logical                        :: inuse
    integer, parameter             :: nmin=10   ! avoid lower numbers which are sometimes reserved
    integer, parameter             :: nmax=999  ! may be system-dependent

    do n = nmin, nmax
       inquire(unit=n, opened=inuse)
       if (.not. inuse) then
          if (present(unit)) unit=n
          return
       end if
    end do
    stop "newunit ERROR: all units are in use"
  end function newunit

end module utils
