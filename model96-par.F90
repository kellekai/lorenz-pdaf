program lorenz96_seq

    use mpi
    implicit none

    !include 'mpif.h'
    ! for distributing along the ranks if q = NG/mpi_size not integer
    integer, parameter                              :: MPI_MIN_BLK = 1

    integer, parameter                              :: NG = 64
    integer, parameter                              :: NT = 1000
    double precision, parameter                     :: F  = 0.2
    double precision, parameter                     :: dt = 0.01
    double precision, allocatable, dimension(:)     :: x
    double precision, allocatable, dimension(:)     :: x_old
    double precision, allocatable, dimension(:)     :: ki
    double precision, allocatable, dimension(:)     :: kj
    integer                                         :: i,j
    integer                                         :: ierr

    integer                                         :: mpi_rank
    integer                                         :: mpi_size
    integer                                         :: mpi_left, mpi_right

    integer                                         :: nl, nlt
    integer                                         :: nl_mod
    integer, allocatable, dimension(:)              :: nl_all

    integer                                         :: dbg_var_int

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

    allocate( x(nlt) )   
    allocate( x_old(nlt) )   
    allocate( ki(nlt) )
    allocate( kj(nlt) )

    x = 0.0d0
    if (mpi_rank .eq. 0) then
        x(3) = 1.0d0
    endif
    do i = 1, nt
        ! use runga kutte RK4 method to solve lorenz96
        ! https://en.wikipedia.org/wiki/Runge-Kutta_methods 
        x_old   = x
        ki      = x
        call d96(ki, kj, F)
        x       = x + dt * kj/6.0
        ki      = x_old + dt * kj/2.0
        if (modulo(mpi_rank,2) .eq. 0) then
            call mpi_send(ki(nlt-2), 2, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, ierr)
            call mpi_recv(ki(1), 2, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, &
                MPI_STATUS_IGNORE, ierr)
            call mpi_send(ki(3), 1, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, ierr)
            call mpi_recv(ki(nlt), 1, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, &
                MPI_STATUS_IGNORE, ierr)
        else
            call mpi_recv(ki(1), 2, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, &
                MPI_STATUS_IGNORE, ierr)
            call mpi_send(ki(nlt-2), 2, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, ierr)
            call mpi_recv(ki(nlt), 1, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, &
                MPI_STATUS_IGNORE, ierr)
            call mpi_send(ki(3), 1, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, ierr)
        endif
        call d96(ki, kj, F)
        x       = x + dt * kj/3.0
        ki      = x_old + dt * kj/2.0 
        if (modulo(mpi_rank,2) .eq. 0) then
            call mpi_send(ki(nlt-2), 2, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, ierr)
            call mpi_recv(ki(1), 2, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, &
                MPI_STATUS_IGNORE, ierr)
            call mpi_send(ki(3), 1, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, ierr)
            call mpi_recv(ki(nlt), 1, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, &
                MPI_STATUS_IGNORE, ierr)
        else
            call mpi_recv(ki(1), 2, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, &
                MPI_STATUS_IGNORE, ierr)
            call mpi_send(ki(nlt-2), 2, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, ierr)
            call mpi_recv(ki(nlt), 1, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, &
                MPI_STATUS_IGNORE, ierr)
            call mpi_send(ki(3), 1, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, ierr)
        endif
        call d96(ki, kj, F)
        x       = x + dt * kj/3.0
        ki      = x_old + dt * kj 
        if (modulo(mpi_rank,2) .eq. 0) then
            call mpi_send(ki(nlt-2), 2, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, ierr)
            call mpi_recv(ki(1), 2, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, &
                MPI_STATUS_IGNORE, ierr)
            call mpi_send(ki(3), 1, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, ierr)
            call mpi_recv(ki(nlt), 1, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, &
                MPI_STATUS_IGNORE, ierr)
        else
            call mpi_recv(ki(1), 2, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, &
                MPI_STATUS_IGNORE, ierr)
            call mpi_send(ki(nlt-2), 2, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, ierr)
            call mpi_recv(ki(nlt), 1, MPI_DOUBLE_PRECISION, mpi_right, 42, MPI_COMM_WORLD, &
                MPI_STATUS_IGNORE, ierr)
            call mpi_send(ki(3), 1, MPI_DOUBLE_PRECISION, mpi_left, 42, MPI_COMM_WORLD, ierr)
        endif
        call d96(ki, kj, F)
        x       = x + dt * kj/6.0 
        ! exchange halo
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
    end do
    

    dbg_var_int = 0
    do i = 1, mpi_rank
        dbg_var_int = dbg_var_int + nl_all(i)
    end do
    do i = 1, nl
        print *, dbg_var_int + i, x(i+2)
    end do
    
    deallocate(nl_all)
    deallocate(x)
    deallocate(x_old)
    deallocate(ki)
    deallocate(kj)

    call mpi_finalize(ierr)

contains
    subroutine d96(x, d, F)
        double precision, dimension(:), intent(IN)      :: x
        double precision, dimension(:), intent(OUT)     :: d
        double precision, intent(IN)                    :: F
        integer                                         :: N
        integer                                         :: i

        N = size(x)
        do i = 3,N-1
            d(i) = ( x(i+1) - x(i-2) ) * x(i-1) - x(i)
        end do
        d = d + F
    end subroutine

end program lorenz96_seq
