program lorenz96_seq
  implicit none
  
  integer, parameter                        :: N    = 30
  integer, parameter                        :: NT   = 1000
  real(8), parameter                        :: F    = 0.2
  real(8), parameter                        :: dt   = 0.01
  real(8), allocatable, dimension(:)        :: x
  real(8), allocatable, dimension(:)        :: d
  integer                                   :: ti
   
  allocate( x(N) )   
  allocate( d(N) )   

  x = 0.0d0
  do ti = 1, NT
      call d96(x, d, F)
      x = x + d * dt
  end do

  do ti = 1, N
      print *, x(ti)
  end do

  contains
    subroutine d96(x, d, F)
        real(8), dimension(:), intent(IN)   :: x
        real(8), dimension(:), intent(OUT)  :: d
        real(8), intent(IN)                 :: F
        integer                             :: N
        integer                             :: i

        N = size(x)
        d(N) = ( x(1) - x(N-2) ) * x(N-1) - x(N)
        d(1) = ( x(2) - x(N-1) ) * x(N)   - x(1) 
        d(2) = ( x(3) - x(N)   ) * x(1)   - x(2)
        do i = 3,N-1
            d(i) = ( x(i+1) - x(i-2) ) * x(i-1) - x(i)
        end do
        d = d + F
    end subroutine

end program lorenz96_seq
