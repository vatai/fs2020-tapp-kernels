!-------------------------------------------------------------------------------
!
!+  Program common kernel driver (streamlike)
!
!-------------------------------------------------------------------------------

#include "profiler.h"

program commonkernel_streamlike
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_streamlike, only: &
     streamlike_pattern1, &
     streamlike_pattern2, &
     streamlike_pattern3
  !-----------------------------------------------------------------------------
  implicit none

  real(RP),dimension(16900,96,1) :: metrics
  real(RP),dimension(16900,96,1,6) :: PROG
  real(RP),dimension(16900,96,1,6) :: PROGq
  real(RP),dimension(16900,96,1) :: rho
  real(RP),dimension(16900,96,1) :: vx
  real(RP),dimension(16900,96,1) :: vy
  real(RP),dimension(16900,96,1) :: vz
  real(RP),dimension(16900,96,1) :: ein
  real(RP),dimension(16900,96,1,6) :: q

  real(RP),dimension(16900,96,1,6) :: tendency_0
  real(RP),dimension(16900,96,1) :: tendency_11
  real(RP),dimension(16900,96,1) :: tendency_12
  real(RP),dimension(16900,96,1) :: tendency_13
  real(RP),dimension(16900,96,1) :: tendency_14
  real(RP),dimension(16900,96,1) :: tendency_15
  real(RP),dimension(16900,96,1) :: tendency_16
  real(RP),dimension(16900,96,1) :: tendency_21
  real(RP),dimension(16900,96,1) :: tendency_22
  real(RP),dimension(16900,96,1) :: tendency_23
  real(RP),dimension(16900,96,1) :: tendency_24
  real(RP),dimension(16900,96,1) :: tendency_25
  real(RP),dimension(16900,96,1) :: tendency_26
  real(RP),dimension(16900,96,1,6) :: tendency

  real(RP),dimension(16900,96,1,6) :: value
  real(RP)              :: fraction

  integer  :: iteration
  !=============================================================================

  metrics(:,:,:)   = 1.0_RP / 3.0_RP
  PROG   (:,:,:,:) = 1.0_RP / 3.0_RP
  PROGq  (:,:,:,:) = 1.0_RP / 3.0_RP

  tendency_0 (:,:,:,:) = 1.0_RP / 3.0_RP
  tendency_11(:,:,:)   = 1.0_RP / 3.0_RP
  tendency_12(:,:,:)   = 1.0_RP / 3.0_RP
  tendency_13(:,:,:)   = 1.0_RP / 3.0_RP
  tendency_14(:,:,:)   = 1.0_RP / 3.0_RP
  tendency_15(:,:,:)   = 1.0_RP / 3.0_RP
  tendency_16(:,:,:)   = 1.0_RP / 3.0_RP
  tendency_21(:,:,:)   = 1.0_RP / 3.0_RP
  tendency_22(:,:,:)   = 1.0_RP / 3.0_RP
  tendency_23(:,:,:)   = 1.0_RP / 3.0_RP
  tendency_24(:,:,:)   = 1.0_RP / 3.0_RP
  tendency_25(:,:,:)   = 1.0_RP / 3.0_RP
  tendency_26(:,:,:)   = 1.0_RP / 3.0_RP

  value(:,:,:,:) = 1.0_RP / 3.0_RP

  !###############################################################################

  PROF_INIT
  PROF_START_ALL
  PROF_START(1)
!  do iteration = 1, SET_iteration ! 100 times
  do iteration = 1, 1 ! 100 times

     call streamlike_pattern2( 16900,                & ! [IN]
                               96,                   & ! [IN]
                               1,                    & ! [IN]
                               tendency_0 (:,:,:,:), & ! [IN]
                               tendency_11(:,:,:),   & ! [IN]
                               tendency_12(:,:,:),   & ! [IN]
                               tendency_13(:,:,:),   & ! [IN]
                               tendency_14(:,:,:),   & ! [IN]
                               tendency_15(:,:,:),   & ! [IN]
                               tendency_16(:,:,:),   & ! [IN]
                               tendency_21(:,:,:),   & ! [IN]
                               tendency_22(:,:,:),   & ! [IN]
                               tendency_23(:,:,:),   & ! [IN]
                               tendency_24(:,:,:),   & ! [IN]
                               tendency_25(:,:,:),   & ! [IN]
                               tendency_26(:,:,:),   & ! [IN]
                               tendency   (:,:,:,:)  ) ! [OUT]

  enddo

  PROF_STOP(1)
  PROF_STOP_ALL
  PROF_FINALIZE

! stop
end program commonkernel_streamlike
!-------------------------------------------------------------------------------
