!$Id$
!BOP
!
! !ROUTINE: obs_op_pdaf --- Implementation of observation operator
!
! !INTERFACE:
SUBROUTINE obs_op_pdaf(step, dim_p, dim_obs_p, state_p, m_state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/ETKF/ESTKF
!
! The routine is called during the analysis step.
! It has to perform the operation of the
! observation operator acting on a state vector.
! For domain decomposition, the action is on the
! PE-local sub-domain of the state and has to 
! provide the observed sub-state for the PE-local 
! domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
	USE mod_assimilation, &
	     ONLY: obs_index_p, obs_p, epoch

  USE mod_parallel, &
    ONLY: mype_filter
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step               ! Currrent time step
  INTEGER, INTENT(in) :: dim_p              ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_p          ! Dimension of observed state
  REAL, INTENT(in)    :: state_p(dim_p)     ! PE-local model state
  REAL, INTENT(out) :: m_state_p(dim_obs_p) ! PE-local observed state

! !CALLING SEQUENCE:
! Called by: PDAF_seek_analysis   (as U_obs_op)
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT
! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
!EOP

! *** local variables ***
	INTEGER :: i
  CHARACTER(len=5)                                :: mpestr
  CHARACTER(len=5)                                :: epostr
! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

  write(mpestr,'(i5.5)') mype_filter
  write(epostr,'(i5.5)') epoch
  open(10, file='../output/obs_server_rank'//TRIM(mpestr)//'_epoch'//TRIM(epostr)//'.txt', form='formatted')
  DO i = 1, dim_obs_p
    m_state_p(i) = state_p(obs_index_p(i))
    write(10,"(I4, 1X, F18.17, 1X, F18.17)") obs_index_p(i), obs_p(i), m_state_p(i)
  END DO
  close(10)

END SUBROUTINE obs_op_pdaf
