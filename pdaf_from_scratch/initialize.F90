!$Id$
!BOP
!
! !ROUTINE: initialize  --- initialize the model for PDAF offline
!
! !INTERFACE:
SUBROUTINE initialize()

! !DESCRIPTION:
! Routine to perform initialization of the model infromation for
! PDAF. Here, the global size of the model domain, the global size
! of the model state vector and the sizes for decomposition of the 
! state vector need to be initialized.
! Generally, this could also be joined with the routine init_pdaf().
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &        ! Model variables
        ONLY: dim_state_p, dim_state, state_min_p, state_max_p
  USE mod_parallel, &     ! Parallelization variables
       ONLY: mype_world, mype_filter, npes_filter, MPI_BLK_DECO

  IMPLICIT NONE

! !ARGUMENTS:
!EOP

! local variables
  INTEGER, ALLOCATABLE  :: dim_state_all(:)
  INTEGER               :: dim_mod
  INTEGER               :: i
! *** Model specifications ***

! nx, ny = ?
  
  ALLOCATE(dim_state_all(npes_filter))

  dim_state = 1024
  dim_state_all = dim_state / npes_filter
  dim_mod = modulo( dim_state, npes_filter )
  do while (dim_mod .gt. 0)
      do i = 1, npes_filter
          if (dim_mod .gt. MPI_BLK_DECO) then
              dim_state_all(i) = dim_state_all(i) + MPI_BLK_DECO
              dim_mod = dim_mod - MPI_BLK_DECO
          else
              dim_state_all(i) = dim_state_all(i) + dim_mod
              dim_mod = 0
              exit
          endif
      end do
  end do
 
  state_min_p = 1
  do i=1,mype_filter
    state_min_p = state_min_p + dim_state_all(i)
  end do
  dim_state_p = dim_state_all(mype_filter+1)
  state_max_p = state_min_p + dim_state_p - 1
  
!  print *, 'rank', mype_filter
!  print *, 'dim_state_p', dim_state_p
!  print *, 'state_min_p', state_min_p
!  print *, 'state_max_p', state_max_p

! *** Screen output ***
  screen2: IF (mype_world == 0) THEN
     WRITE (*, '(1x, a)') 'INITIALIZE MODEL INFORMATION FOR PDAF OFFLINE MODE'
     WRITE (*, '(22x,a)') 'MODEL: ...'
!      WRITE (*, '(24x,a,i4,1x,a1,1x,i4)') 'Grid size:',nx,'x',ny
     WRITE (*, '(5x, a, i7)') &
          'Global model state dimension:', dim_state
  END IF screen2

! *** Initialize dimensions and fields with domain decomposition



END SUBROUTINE initialize
