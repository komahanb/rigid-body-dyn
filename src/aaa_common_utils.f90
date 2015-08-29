module common_utils
  use constants
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

  !---------------------------------------------------------------------------
  ! function that reverses the order of the constituent elements in the vector
  !---------------------------------------------------------------------------
  function reverse_real(a) result(rev_a)
    real(dp), intent(in)     :: a(:)
    real(dp)                 :: rev_a(size(a))
    integer(sp)              :: i , n

    n= size(a)
    do i =  1, n
       rev_a(i) = a(n+1-i)
    end do
    
  end function reverse_real

  !---------------------------------------------------------------------------
  ! function that reverses the order of the constituent elements in the vector
  !---------------------------------------------------------------------------
  function reverse_int(a) result(rev_a)
    integer(sp), intent(in)     :: a(:)
    integer(sp)                 :: rev_a(size(a))
    integer(sp)                 :: i , n

    n= size(a)
    do i =  1, n
       rev_a(i) = a(n+1-i)
    end do
    
  end function reverse_int
  
end module common_utils
