!$Id$
!BOP
!
! !ROUTINE: init_dim_obs_pdaf --- Compute number of observations
!
! !INTERFACE:
SUBROUTINE init_dim_obs_pdaf(step, dim_obs_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF
!
! The routine is called at the beginning of each
! analysis step.  It has to initialize the size of 
! the observation vector according to the current 
! time step for the PE-local domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &
      ONLY :  dim_state, &      ! dimension of full state
              obs_blk_size, &   ! block size for observations
              obs_prcnt, &      ! percentage of full state = number of observations
              state_min_p, &    ! minimum index state on PE
              state_max_p, &    ! maximum index state on PE
              obs_index_p       ! index array of observations
  !USE mod_parallel, &
  !    ONLY: mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step       ! Current time step
  INTEGER, INTENT(out) :: dim_obs_p  ! Dimension of observation vector

! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis    (as U_init_dim_obs)
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
! Called by: PDAF_etkf_analysis, PDAF_etkf_analysis_T
! Called by: PDAF_estkf_analysis, PDAF_estkf_analysis_fixed
! Called by: PDAF_netf_analysis
!EOP

! *** Local variables
  INTEGER               :: i, j                     ! Counters
  INTEGER               :: num_reg    ! number of regions
  INTEGER               :: dim_obs    ! number of observations global
  INTEGER               :: stride     ! stride for the regions
  INTEGER               :: offset     ! offset for index in loop
  INTEGER               :: index_tmp  ! index dummy var
  INTEGER               :: cnt_obs_p  ! counter for observations on PE
! ****************************************
! *** Initialize observation dimension ***
! ****************************************

! compute total number of observations
  dim_obs = obs_prcnt * dim_state
  if (dim_obs .eq. 0) then
    dim_obs = 1
  end if

! compute number of regions
  num_reg = dim_obs / obs_blk_size
  if ( MODULO(dim_obs, obs_blk_size) .ne. 0 ) then
    num_reg = num_reg + 1
  end if

! compute stride for regions
  stride = dim_state / num_reg
  
! determine number of obs in pe
  dim_obs_p = 0
  do i=1,num_reg
    offset = (i-1) * stride
    do j=1,obs_blk_size
      index_tmp = offset + j
      if ( (index_tmp .ge. state_min_p) .and. (index_tmp .le. state_max_p) ) then
        dim_obs_p = dim_obs_p + 1
      end if
      if ( index_tmp .eq. state_max_p ) exit
    end do
    if ( index_tmp .eq. state_max_p ) exit
  end do

  ALLOCATE( obs_index_p(dim_obs_p) )
  
! assign indices to index array
  cnt_obs_p = 0
  do i=1,num_reg
    offset = (i-1) * stride
    do j=1,obs_blk_size
      index_tmp = offset + j
      if ( (index_tmp .ge. state_min_p) .and. (index_tmp .le. state_max_p) ) then
        cnt_obs_p = cnt_obs_p + 1
        obs_index_p(cnt_obs_p) = index_tmp
      end if
      if ( index_tmp .eq. state_max_p ) exit
    end do
    if ( index_tmp .eq. state_max_p ) exit
  end do


  ! Template reminder - delete when implementing functionality
  WRITE (*,*) 'TEMPLATE init_dim_obs_pdaf.F90: Initialize observation dimension here!'


! dim_obs_p = ?

END SUBROUTINE init_dim_obs_pdaf

