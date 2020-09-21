program lorenz96_seq

  use mpi
  USE parser, &           ! Parser function
    ONLY: parse
  implicit none

  !include 'mpif.h'
  ! for distributing along the ranks if q = NG/mpi_size not integer
  integer, parameter                              :: MPI_MIN_BLK = 1

  integer, parameter                              :: NG = 1024
  integer, parameter                              :: NT = 1000
  real, parameter                     :: F  = 0.2
  real, parameter                     :: dt = 0.01
  real, allocatable, dimension(:)     :: x
  real, allocatable, dimension(:)     :: x_old
  real, allocatable, dimension(:)     :: ki
  real, allocatable, dimension(:)     :: kj
  integer                                         :: i,j
  integer                                         :: ierr

  integer                                         :: mpi_rank
  integer                                         :: mpi_size
  integer                                         :: mpi_left, mpi_right
  
  integer                                         :: nl, nlt
  integer                                         :: nl_mod
  integer, allocatable, dimension(:)              :: nl_all
  integer(8)                                      :: nl_off
  integer                                         :: state_min_p
  integer                                         :: state_max_p

  integer                                         :: obs_block = 32
  integer                                         :: obs_share = 5
  real                                            :: obs_percent
  integer                                         :: dbg_var_int
  integer, parameter                              :: rnd_blk_size =32
  integer                                         :: member = 1
  CHARACTER(len=5)                                :: ensstr
  CHARACTER(len=5)                                :: memstr
  CHARACTER(len=5)                                :: epostr
  integer(4)                                      :: epoch = 0
  real           :: mean = 0.0d0
  real           :: stdv = 0.2d0
  integer(4), parameter :: seed_init = 310780
  integer(4)        :: seed
  CHARACTER(len=32) :: handle  ! handle for command line parser
  integer :: generate_obs = 0
  CHARACTER(len=256) :: data_path = ''

  call init_parallel()

  !   parse commandline args
  handle = 'member'             ! Control application of model error
  CALL parse(handle, member)
  handle = 'obs_gen'             ! Control application of model error
  CALL parse(handle, generate_obs)
  handle = 'seed'             ! Control application of model error
  CALL parse(handle, seed)
  handle = 'epoch'             ! Control application of model error
  CALL parse(handle, epoch)
  handle = 'obs_block'             ! Control application of model error
  CALL parse(handle, obs_block)
  handle = 'obs_share'             ! Control application of model error
  CALL parse(handle, obs_share)
  obs_percent = real(obs_share)/100
  !handle = 'data_path'             ! Control application of model error
  !CALL parse(handle, data_path)
  
  data_path = 'pdaf/'

  allocate( x(nlt) )   
  allocate( x_old(nlt) )   
  allocate( ki(nlt) )
  allocate( kj(nlt) )

! create initial state with some noise and a 
! small perturbation
  if ( epoch .eq. 0 ) then
    x = 0.0
    if (mpi_rank .eq. 0) then
      x(5) = 1.0
    end if
    seed = seed_init
    call add_noise( x(3:3+nl-1), nl, 0.02, seed )
  else
    WRITE (ensstr, '(i5.5)') member
    call read_ens(TRIM(data_path)//'ens_'//TRIM(ensstr)//'_ana.txt')
  end if
  
  if (member .eq. 1) then
      write(epostr,'(i5.5)') epoch
      write(memstr,'(i5.5)') member
      open(10, file='output/ana_'//TRIM(epostr)//'_'//TRIM(memstr)//'.txt', form='formatted')
      write(10,"(F18.17)") x(3:3+nl-1)
      close(10)
  end if

  do i = 1, nt
    ! use runga kutte RK4 method to solve lorenz96
    ! https://en.wikipedia.org/wiki/Runge-Kutta_methods 
    x_old   = x
    ki      = x
    call d96(ki, kj, F)
    x       = x + dt * kj/6.0
    ki      = x_old + dt * kj/2.0
    call exchange(ki)
    call d96(ki, kj, F)
    x       = x + dt * kj/3.0
    ki      = x_old + dt * kj/2.0 
    call exchange(ki)
    call d96(ki, kj, F)
    x       = x + dt * kj/3.0
    ki      = x_old + dt * kj 
    call exchange(ki)
    call d96(ki, kj, F)
    x       = x + dt * kj/6.0 
    call exchange(x)
  end do
  
  if (generate_obs .eq. 0) then
      if ( mpi_rank .eq. 0 ) then
        write(epostr,'(i5.5)') epoch
        write(memstr,'(i5.5)') member
        open(10, file='output/for_'//TRIM(epostr)//'_'//TRIM(memstr)//'.txt', form='formatted')
        write(10,"(F18.17)") x(3:3+nl-1)
        close(10)
      end if
    WRITE (ensstr, '(i5.5)') member
    call write_ens(TRIM(data_path)//'ens_'//TRIM(ensstr)//'_for.txt')
  else
!   simulate a small error in time direction
!   only if number of timesteps is large enough
    do i = 1,int(0.01*nt)
      ! use runga kutte RK4 method to solve lorenz96
      ! https://en.wikipedia.org/wiki/Runge-Kutta_methods 
      x_old   = x
      ki      = x
      call d96(ki, kj, F)
      x       = x + dt * kj/6.0
      ki      = x_old + dt * kj/2.0
      call exchange(ki)
      call d96(ki, kj, F)
      x       = x + dt * kj/3.0
      ki      = x_old + dt * kj/2.0 
      call exchange(ki)
      call d96(ki, kj, F)
      x       = x + dt * kj/3.0
      ki      = x_old + dt * kj 
      call exchange(ki)
      call d96(ki, kj, F)
      x       = x + dt * kj/6.0 
      call exchange(x)
    end do
    
    seed = seed_init
    call add_noise( x(3:3+nl-1), nl, 0.01, seed )
    
    call write_obs(TRIM(data_path)//'obs.txt', x(3:3+nl-1), obs_percent, obs_block)
  end if

  deallocate(nl_all)
  deallocate(x)
  deallocate(x_old)
  deallocate(ki)
  deallocate(kj)

  500 call mpi_finalize(ierr)

contains
  subroutine d96(x, d, F)
    real, dimension(:), intent(IN)      :: x
    real, dimension(:), intent(OUT)     :: d
    real, intent(IN)                    :: F
    integer                                         :: N
    integer                                         :: i

    N = size(x)
    do i = 3,N-1
      d(i) = ( x(i+1) - x(i-2) ) * x(i-1) - x(i)
    end do
    d = d + F
  end subroutine

  subroutine add_noise( x, n, sigma, seed )
    real, intent(inout), dimension(:)   :: x
    integer, intent(in)                             :: n
    real, intent(in)                             :: sigma
    integer, intent(in)                             :: seed
    real                                         :: mean = 0.0
    integer(4), parameter                           :: buf_size = 1024
    real, dimension(buf_size)                    :: buf
    integer                                         :: i,j,k

    do i=0,mpi_size-1
      if ( i .eq. mpi_rank ) then
        if ( mpi_rank .gt. 0 ) then
          call mpi_recv( seed, 1, MPI_INTEGER, mpi_rank-1, 0, &
            mpi_comm_world, mpi_status_ignore, ierr ) 
        endif

        do j=1,n
          if(modulo(j, buf_size) .eq. 1) then
            call r8vec_normal_ab(buf_size, mean, sigma, seed, buf )
          endif
          k = modulo(j-1, buf_size) + 1
          x(j) = x(j) + buf(k)
        end do

        if ( i .lt. (mpi_size-1) ) then
          call mpi_send( seed, 1, MPI_INTEGER, mpi_rank+1, 0, &
            mpi_comm_world, ierr ) 
        endif
      end if
    end do


  end subroutine

  subroutine init_parallel()
    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, mpi_size, ierr)

    allocate( nl_all(mpi_size) )

    ! middle ranks
    mpi_left    = mpi_rank - 1
    mpi_right   = mpi_rank + 1

    ! first and last rank
    if (mpi_rank == 0) then 
      mpi_left = mpi_size - 1
    elseif (mpi_rank == mpi_size-1) then
      mpi_right = 0
    endif

    nl_all = NG / mpi_size
    nl_mod = modulo( NG, mpi_size )
    do while (nl_mod .gt. 0)
      do i = 1, mpi_size
        if (nl_mod .gt. MPI_MIN_BLK) then
          nl_all(i) = nl_all(i) + MPI_MIN_BLK
          nl_mod = nl_mod - MPI_MIN_BLK
        else
          nl_all(i) = nl_all(i) + nl_mod
          nl_mod = 0
          exit
        endif
      end do
    end do

    nl = nl_all(mpi_rank+1)
    nlt = nl + 3
    
    nl_off = 0
    do i=1,mpi_rank
      nl_off = nl_off + nl_all(i)
    end do
    state_min_p = nl_off + 1
    state_max_p = nl_off + nl

  end subroutine

  subroutine exchange(x)
    real, dimension(:), intent(INOUT)      :: x

    if (modulo(mpi_rank,2) .eq. 0) then
      call mpi_send(x(nlt-2), 2, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, ierr)
      call mpi_recv(x(1), 2, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, &
        MPI_STATUS_IGNORE, ierr)
      call mpi_send(x(3), 1, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, ierr)
      call mpi_recv(x(nlt), 1, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, &
        MPI_STATUS_IGNORE, ierr)
    else
      call mpi_recv(x(1), 2, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, &
        MPI_STATUS_IGNORE, ierr)
      call mpi_send(x(nlt-2), 2, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, ierr)
      call mpi_recv(x(nlt), 1, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, &
        MPI_STATUS_IGNORE, ierr)
      call mpi_send(x(3), 1, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, ierr)
    endif
  end subroutine

  subroutine write_ens( file_name )
    integer     :: thefile
    integer(kind=MPI_OFFSET_KIND) :: disp
    integer     :: ierr
    CHARACTER(len=*) :: file_name

    call mpi_file_open(MPI_COMM_WORLD, file_name, & 
      MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
      MPI_INFO_NULL, thefile, ierr)  

    disp = 0

    do i = 1, mpi_rank
      disp = disp + nl_all(i) * sizeof( x(3) )
    end do

    call mpi_file_set_view(thefile, disp, MPI_DOUBLE_PRECISION, & 
      MPI_DOUBLE_PRECISION, 'native', & 
      MPI_INFO_NULL, ierr)

    call mpi_file_write(thefile, x(3), nl, MPI_DOUBLE_PRECISION, & 
      MPI_STATUS_IGNORE, ierr)

    call mpi_file_close(thefile, ierr)
  end subroutine
  
  subroutine read_ens( file_name )
    integer     :: thefile
    integer(kind=MPI_OFFSET_KIND) :: disp
    integer     :: ierr
    CHARACTER(len=*) :: file_name

    call mpi_file_open(MPI_COMM_WORLD, file_name, & 
      MPI_MODE_RDONLY, & 
      MPI_INFO_NULL, thefile, ierr)  

    disp = 0

    do i = 1, mpi_rank
      disp = disp + nl_all(i) * sizeof( x(3) )
    end do

    call mpi_file_set_view(thefile, disp, MPI_DOUBLE_PRECISION, & 
      MPI_DOUBLE_PRECISION, 'native', & 
      MPI_INFO_NULL, ierr)

    call mpi_file_read(thefile, x(3), nl, MPI_DOUBLE_PRECISION, & 
      MPI_STATUS_IGNORE, ierr)

    call mpi_file_close(thefile, ierr)
  end subroutine

  subroutine write_obs( file_name, x, share, blk_size )
    CHARACTER(len=*)  :: file_name
    real, dimension(:), intent(in) :: x
    real, intent(in)                :: share
    integer, intent(in)             :: blk_size
    integer                         :: dim_obs
    integer                         :: num_reg
    integer                         :: stride
    integer                         :: dim_obs_p
    real, allocatable, dimension(:) :: obs_p
    integer                         :: offset
    integer                         :: index_tmp
    integer                         :: cnt_obs_p
    integer                         :: cnt_obs
    integer                         :: dim_obs_all(mpi_size)
    
    integer     :: thefile
    integer(kind=MPI_OFFSET_KIND) :: disp

    ! compute total number of observations
    dim_obs = share * NG
    if (dim_obs .eq. 0) then
      dim_obs = 1
    end if

    ! compute number of regions
    num_reg = dim_obs / blk_size
    if ( MODULO(dim_obs, blk_size) .ne. 0 ) then
      num_reg = num_reg + 1
    end if

    ! compute stride for regions
    stride = NG / num_reg

    ! determine number of obs in pe
    dim_obs_p = 0
    cnt_obs = 0
    do i=1,num_reg
      offset = (i-1) * stride
      do j=1,blk_size
        index_tmp = offset + j
        if ( (index_tmp .ge. state_min_p) .and. (index_tmp .le. state_max_p) ) then
          dim_obs_p = dim_obs_p + 1
        end if
        cnt_obs = cnt_obs + 1
        if ( cnt_obs .eq. dim_obs ) exit
        if ( index_tmp .eq. state_max_p ) exit
      end do
      if ( cnt_obs .eq. dim_obs ) exit
      if ( index_tmp .eq. state_max_p ) exit
    end do

    ALLOCATE( obs_p(dim_obs_p) )

    ! assign indices to index array
    cnt_obs_p = 0
    do i=1,num_reg
      offset = (i-1) * stride
      do j=1,blk_size
        index_tmp = offset + j
        if ( (index_tmp .ge. state_min_p) .and. (index_tmp .le. state_max_p) ) then
          cnt_obs_p = cnt_obs_p + 1
          obs_p(cnt_obs_p) = x(index_tmp - (state_min_p - 1))
        end if
        if ( cnt_obs_p .eq. dim_obs_p ) exit
      end do
      if ( cnt_obs_p .eq. dim_obs_p ) exit
    end do

    ! write observations
    call mpi_allgather( dim_obs_p, 1, MPI_INTEGER, dim_obs_all, &
      1, MPI_INTEGER, mpi_comm_world, ierr )

    disp = 0

    do i = 1, mpi_rank
      disp = disp + dim_obs_all(i) * sizeof( obs_p(1) )
    end do

    call mpi_file_open(MPI_COMM_WORLD, file_name, & 
      MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
      MPI_INFO_NULL, thefile, ierr)  

    call mpi_file_set_view(thefile, disp, MPI_DOUBLE_PRECISION, & 
      MPI_DOUBLE_PRECISION, 'native', & 
      MPI_INFO_NULL, ierr)

    call mpi_file_write(thefile, obs_p, dim_obs_p, MPI_DOUBLE_PRECISION, & 
      MPI_STATUS_IGNORE, ierr)

    call mpi_file_close(thefile, ierr)
    

    if ( mpi_rank .eq. 0 ) then
      write(epostr,'(i5.5)') epoch
      open(10, file='output/obs_'//TRIM(epostr)//'.txt', form='formatted')
      write(10,"(F18.17)") obs_p
      close(10)
    end if


  end subroutine

end program lorenz96_seq
