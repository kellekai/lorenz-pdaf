!$Id$
!BOP
!
! !ROUTINE: init_ens_offline --- Initialize ensemble for SEIK in offline mode
!
! !INTERFACE:
SUBROUTINE init_ens_offline(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! The routine is called when the filter is
! initialized in PDAF\_filter\_init.  It has
! to initialize an ensemble of dim\_ens states.
! For the offline mode, the ensemble will be
! typically read-in from files.
!
! The routine is called by all filter processes and 
! initializes the ensemble for the PE-local domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code based on offline_1D
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
       ONLY: state_min_p
  USE mod_parallel, &
       ONLY: mype_filter, COMM_filter, MPI_DOUBLE_PRECISION, &
       MPI_INFO_NULL, MPI_MODE_RDONLY, MPI_STATUS_IGNORE

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state
  ! It is not necessary to initialize the array 'state_p' for SEIK. 
  ! It is available here only for convenience and can be used freely.
  REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag
  INTEGER       :: file_id      ! MPI file handle
  INTEGER       :: ierr         ! MPI error handle

! !CALLING SEQUENCE:
! Called by: PDAF_filter_init    (as U_ens_init)
!EOP

! *** local variables ***
   INTEGER :: i, j, col, member            ! Counters
   INTEGER, SAVE :: allocflag = 0     ! Flag for memory counting
   REAL, ALLOCATABLE :: ens(:,:)      ! global ensemble array
   REAL, ALLOCATABLE :: field(:,:)    ! global model field
   CHARACTER(len=2) :: ensstr         ! String for ensemble member
   ! variables and arrays for domain decomposition
   INTEGER :: offset                  ! Row-offset according to domain decomposition
   INTEGER :: domain                  ! domain counter
   REAL,ALLOCATABLE :: ens_p_tmp(:,:) ! Temporary ensemble for some PE-domain


! **********************
! *** INITIALIZATION ***
! **********************

  mype0a: IF (mype_filter==0) THEN

     WRITE (*, '(/9x, a)') 'Initialize state ensemble'
     WRITE (*, '(9x, a)') '--- read ensemble from files'
     WRITE (*, '(9x, a, i5)') '--- Ensemble size:  ', dim_ens

  END IF mype0a


! ********************************
! *** Read ensemble from files ***
! ********************************

  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE init_ens_offline.F90: Initialize ensemble array here!'
    
  do member=1, dim_ens
		WRITE (ensstr, '(i3.3)') member
    call MPI_FILE_OPEN(COMM_filter, 'ens_'//TRIM(ensstr)//'_for.txt', & 
                     MPI_MODE_RDONLY, & 
                     MPI_INFO_NULL, file_id, ierr) 
    call MPI_FILE_SET_VIEW(file_id, state_min_p*sizeof(1.0) , MPI_DOUBLE_PRECISION, & 
                         MPI_DOUBLE_PRECISION, 'native', & 
                         MPI_INFO_NULL, ierr) 
    call MPI_FILE_READ(file_id, ens_p(1,member), dim_p, MPI_DOUBLE_PRECISION, & 
                      MPI_STATUS_IGNORE, ierr) 
    call MPI_FILE_CLOSE(file_id, ierr)   
  end do

  ! ens_p = ?


! ****************************
! *** Distribute substates ***
! ****************************

  ! This is an example how one could distribute ensemble information over multiple processes


!   mype0c: IF (mype_filter == 0) THEN
!      ! *** Initialize and send sub-state on PE 0 ***
! 
!      ! Initialize sub-ensemble for PE 0
!      DO col = 1, dim_ens
!         DO i=1, dim_p
!            ens_p(i, col) = ens(i, col)
!         END DO
!      END DO
! 
!      ! Define offset in state vectors
!      offset = local_dims(1)
! 
!      DO domain = 2, npes_filter
!         ! Initialize sub-ensemble for other PEs and send sub-arrays
! 
!         ! Allocate temporary buffer array
!         ALLOCATE(ens_p_tmp(local_dims(domain), dim_ens))
! 
!         ! Initialize MPI buffer for local ensemble
!         DO col = 1, dim_ens
!            DO i = 1, local_dims(domain)
!               ens_p_tmp(i, col) = ens(i + offset, col)
!            END DO
!         END DO
! 
!         ! Send sub-arrays
!         CALL MPI_send(ens_p_tmp, dim_ens * local_dims(domain), &
!              MPI_DOUBLE_PRECISION, domain - 1, 1, COMM_filter, MPIerr)
! 
!         DEALLOCATE(ens_p_tmp)
! 
!         ! Increment offset
!         offset = offset + local_dims(domain)
! 
!      END DO
! 
!   ELSE mype0c
!      ! *** Receive ensemble substates on filter-PEs with rank > 0 ***
! 
!      CALL MPI_recv(ens_p, dim_p * dim_ens, MPI_DOUBLE_PRECISION, &
!           0, 1, COMM_filter, MPIstatus, MPIerr)
!      
!   END IF mype0c


! ****************
! *** clean up ***
! ****************

!   IF (mype_filter==0) DEALLOCATE(field, ens)

END SUBROUTINE init_ens_offline
