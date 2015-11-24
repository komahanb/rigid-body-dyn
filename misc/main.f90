! program to find the moment of inertia of 3d objects numerically`
program main

  use mpi

  implicit none

  integer, parameter :: nsamples=100000, ndim=3 ! geometry
  real(8)  :: r(nsamples, ndim) ! samples to fill the domain
  real(8)  :: Iglb(ndim) ! mass moment of inertia
  integer  :: j, k, kk ! loop counters
  real(8)  :: I(ndim) ! mass moment of inertia
  real(8)  :: O(ndim) ! location of the axis of rotation

  ! mpi variables
  integer  :: ierr, nproc, idproc
  logical  :: master= .false.
  integer  :: idec, is, ie, id,n, ii, jj
  real(8)  :: bound(ndim,2)
  !*******************************************************************!

  ! MPI initialize
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,idproc,ierr)
  call mpi_comm_size(mpi_comm_world,nproc,ierr)

  ! Identify master
  if (idproc .eq. 0) master = .true.

  ! Generate random numbers in ndim dimensions between 0 to 1
  if (master) then
     call random_seed()
     call random_number(r)
  end if

  ! Broad cast the random number to all procs

  call MPI_BCAST(R,nsamples*ndim,MPI_DOUBLE_PRECISION,0,&
       &MPI_COMM_WORLD,ierr)

  ! Do scaling if needed according to the geometry

  if (master) then
     write(*,'(a)')"==============================================="
     write(*,'(a)')"       Mass Moment of Inertia Calculator       "
     write(*,'(a)')"==============================================="
  end if

  !===============     Main Monte Carlo loop  ========================!

  idec = int(dble(nsamples)/dble(nproc))
  is   = idec*idproc + 1
  ie   = idec*(idproc+1)
  if(idproc .eq.  nproc-1) ie = nsamples ! last proc

  write(*,'(11x,a,2i10,a,i3)') '>> [',is,ie,' ] for Processor',idproc

  ! print*, R, idproc
  bound(1,1) = 0.0
  bound(1,2) = 1.0

  bound(2,1) = 0.0
  bound(2,2) = 0.5

  bound(3,1) = 0.0
  bound(3,2) = 0.1

  ! Scale it to the domain boundaries
  do ii = is, ie
     do jj = 1,ndim
        R(ii,jj) = bound(jj,1)+(bound(jj,2)-bound(jj,1))*R(ii,jj)
     end do
  end do

!!$    
!!$     print*, maxval(R(:,1)), minval(R(:,1)) ! x max, x min
!!$     print*, maxval(R(:,2)), minval(R(:,2)) ! y max, y min
!!$     print*, maxval(R(:,3)), minval(R(:,3)) ! z max, z min

  ! local copy of I to each proc
  O = 0.0d0 ! location of the axis of rotation
  I = 0.0d0

  samples: do k = is, ie
     axis: do kk = 1 , ndim
        dimension: do j = 1, ndim
           if (j .ne. kk)  I(j) = I(j) + (O(j) - r(k,j))**2 
        end do dimension
     end do axis
  end do samples

!!$  ! information sharing
!!$  do id = 0, nproc - 1
!!$     is = idec*id + 1
!!$     ie = idec*(id+1)
!!$     if(id .eq. nproc-1) ie = nsamples
!!$     ! Exchange response values for nsamples
!!$     !     call MPI_BCAST(R, ie-is+1, MPI_DOUBLE_PRECISION, id,&
!!$     !          & MPI_COMM_WORLD, ierr)
!!$  end do

  call mpi_barrier(MPI_COMM_WORLD, IERR)
  call MPI_ALLREDUCE(I,Iglb,ndim,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,ierr)

  Iglb = Iglb / dble(nsamples)

  if (master) then
     write(*,*) "Moment of inertia: ", Iglb 
  end if


  call mpi_finalize(ierr)

end program main
