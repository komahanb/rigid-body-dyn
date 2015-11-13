!============================================================!
! Module that contains functions and routine related to file
! operations such as input / output/ format conversions etc
!============================================================!
module filehandler

  implicit none

contains

  !***********************************************************!
  ! a function that returns a new unit number for file handling
  !***********************************************************!
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

    stop "ERROR: all file units are in use"

  end function newunit

  !***************************************************************!
  ! Routine that coverts integer value to a string value
  !***************************************************************!
  subroutine i_to_s(intval,string)

    implicit none

    integer           :: idig,intval,ipos,ival,i
    character (len=*) :: string

    string(:) = ' '

    ival = intval

    !  Working from right to left, strip off the digits of the integer
    !  and place them into STRING(1:len ( string )).
    ipos = len(string)

    do while ( ival /= 0 )

       idig = mod( ival, 10 )
       ival = ival / 10
       string(ipos:ipos) = char(idig + 48 )
       ipos = ipos - 1

    end do
    
    !Fill the empties with zeroes.
    do i = 1, ipos
       string(i:i) = '0'
    end do

    return

  end subroutine i_to_s


end module filehandler
