#include "profiler.h"

module mod_adm
  implicit none

  integer, public  :: ADM_lall=1
  integer, public  :: ADM_lall_pl=2
  integer, public  :: ADM_gall=16900
  integer, public  :: ADM_gall_pl=6
  integer, public  :: ADM_kall=96
  integer, public  :: ADM_gall_1d=130
#ifdef SINGLE_SIM
  integer, public  :: ADM_gall_1d_in=65
#else
  integer, public  :: ADM_gall_1d_in=130
#endif
  integer, public  :: ADM_gmin=2
  integer, public  :: ADM_gmax=129
  integer, public  :: ADM_gmax_pl=6
  integer, public  :: ADM_kmin=2
  integer, public  :: ADM_kmax=95
  logical, public  :: ADM_have_sgp(1)=.true.
  integer, public,parameter  :: ADM_TI = 1
  integer, public,parameter  :: ADM_TJ = 2
  integer, public,parameter  :: ADM_AI  = 1
  integer, public,parameter  :: ADM_AIJ = 2
  integer, public,parameter  :: ADM_AJ  = 3
  integer, public,parameter  :: ADM_lall_para=1
  integer, public,parameter  :: ADM_gall_para=16900
  integer, public,parameter  :: ADM_kall_para=96
  integer, public,parameter  :: ADM_gall_1d_para=130
#ifdef SINGLE_SIM
  integer, public,parameter  :: ADM_gall_1d_para_in=65
#else
  integer, public,parameter  :: ADM_gall_1d_para_in=130
#endif
end module mod_adm

module mod_precision
  implicit none
  private
  integer, public, parameter :: SP = kind(0.E0) ! Single Precision: kind(0.E0)
  integer, public, parameter :: DP = kind(0.D0) ! Double Precision: kind(0.D0)

#ifdef SINGLE
  integer, public, parameter :: RP = SP
#else
  integer, public, parameter :: RP = DP
#endif
end module mod_precision

module mod_oprt
  use mod_precision
  use mod_adm, only: &
    ADM_lall,       &
    ADM_lall_pl,    &
    ADM_gall,       &
    ADM_gall_pl,    &
    ADM_kall,       &
    ADM_gall_1d,    &
    ADM_gall_1d_in,    &
    ADM_gmin,       &
    ADM_gmax,       &
    ADM_kmin,       &
    ADM_kmax,       &
    ADM_have_sgp,   &
    TI  => ADM_TI,  &
    TJ  => ADM_TJ,  &
    AI  => ADM_AI,  &
    AIJ => ADM_AIJ, &
    AJ  => ADM_AJ,  &
    ADM_lall_para,  &
    ADM_gall_para,  &
    ADM_gall_1d_para,  &
    ADM_gall_1d_para_in, &
    ADM_gmax_pl    

  implicit none

  !++ Public parameters & variables
  integer, public, save :: OPRT_nstart
  integer, public, save :: OPRT_nend

  ! < for diffusion operator >
!  real(RP), public, allocatable, save :: cinterp_TN (:,:,:,:)
!  real(RP), public, allocatable, save :: cinterp_HN (:,:,:,:)
!  real(RP), public, allocatable, save :: cinterp_TRA(:,:,:)
!  real(RP), public, allocatable, save :: cinterp_PRA(:,:)
  real(RP),public,save :: cinterp_TN (ADM_gall_1d_para_in,ADM_gall_1d_para,ADM_lall_para,AI:AJ,1:3)
  real(RP),public,save :: cinterp_HN (ADM_gall_1d_para_in,ADM_gall_1d_para,ADM_lall_para,AI:AJ,1:3)
  real(RP),public,save :: cinterp_TRA(ADM_gall_1d_para_in,ADM_gall_1d_para,ADM_lall_para,TI:TJ)
  real(RP),public,save :: cinterp_PRA(ADM_gall_1d_para_in,ADM_gall_1d_para,ADM_lall_para)

contains

  subroutine OPRT_snap_read ( scl, scl_pl, kh, kh_pl, mfact )
    implicit none
    integer :: i

    real(RP), intent(out)  :: scl    (ADM_gall_1d_in,ADM_gall_1d ,ADM_kall,ADM_lall   )
    real(RP), intent(out)  :: scl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out)  :: kh     (ADM_gall_1d_in,ADM_gall_1d   ,ADM_kall,ADM_lall   )
    real(RP), intent(out)  :: kh_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out)  :: mfact
    integer j,k,l,g,gi,gj
    integer suf

    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
!    allocate( cinterp_TN(ADM_gall,ADM_lall,AI:AJ,1:3) )
!    allocate( cinterp_HN(ADM_gall,ADM_lall,AI:AJ,1:3) )
!    allocate( cinterp_TRA(ADM_gall,ADM_lall,TI:TJ) )
!    allocate( cinterp_PRA(      ADM_gall,ADM_lall) )

#ifndef OUTPUTCHECK
        cinterp_TN(1:ADM_gall_1d_in,1:ADM_gall_1d,1:ADM_lall,AI:AJ,1:3) = 1.0_RP
        cinterp_HN(1:ADM_gall_1d_in,1:ADM_gall_1d,1:ADM_lall,AI:AJ,1:3) = 1.0_RP
        cinterp_TRA(1:ADM_gall_1d_in,1:ADM_gall_1d,1:ADM_lall,TI:TJ) = 1.0_RP
        cinterp_PRA(      1:ADM_gall_1d_in,1:ADM_gall_1d,1:ADM_lall) = 1.0_RP

        scl    (1:ADM_gall_1d_in,1:ADM_gall_1d   ,1:ADM_kall,1:ADM_lall   ) = 1.0_RP
        scl_pl (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0_RP
        kh     (1:ADM_gall_1d_in,1:ADM_gall_1d   ,1:ADM_kall,1:ADM_lall   ) = 1.0_RP
        kh_pl  (1:ADM_gall_pl,1:ADM_kall,1:ADM_lall_pl) = 1.0_RP
#else

        do i=1,3
         do j=AI,AJ
          do l=1,ADM_lall
           do gj = 1,ADM_gall_1d
           do gi = 1,ADM_gall_1d_in
             cinterp_TN(gi,gj,l,j,i) = 1.0 / suf(gi,gj) / l / (i * j)
             cinterp_HN(gi,gj,l,j,i) = 2.0 / suf(gi,gj) / l / (i * j)
           enddo
           enddo
          enddo
         enddo
        enddo
        do i=TI,TJ
         do l=1,ADM_lall
          do gj = 1, ADM_gall_1d
          do gi = 1, ADM_gall_1d_in
            cinterp_TRA(gi,gj,l,i) = 3.0 / suf(gi,gj) / l / i
          enddo
          enddo
         enddo
        enddo

       do l = 1 , ADM_lall
         do gj = 1,ADM_gall_1d
         do gi = 1,ADM_gall_1d_in
          cinterp_PRA(gi,gj,l) = 4.0 / suf(gi,gj) / l
         enddo
         enddo
       enddo
       do l = 1 , ADM_lall
       do k = 1 , ADM_kall
         do gj = 1,ADM_gall_1d
         do gi = 1,ADM_gall_1d_in
            g = suf(gi,gj)
            scl(gi,gj,k,l) = 1.1 / g / l /  k
            kh(gi,gj,k,l) = 2.2 / g / l / k
         enddo
         enddo
       enddo
       enddo

#endif
        mfact = 1.0_RP
        return
  end subroutine OPRT_snap_read

  subroutine OPRT_diffusion( &
       dscl, dscl_pl, &
       scl,  scl_pl,  &
       kh,   kh_pl,   &
       mfact          )
!    use mod_gmtr, only: &
!       GMTR_P_var_pl, &
!       GMTR_T_var_pl, &
!       GMTR_A_var_pl
!    use mod_debug, only: DEBUG_rapstart, DEBUG_rapend
    implicit none

!    real(RP), intent(out) :: dscl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: dscl   (ADM_gall_1d_in,ADM_gall_1d ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: dscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
!    real(RP), intent(in)  :: scl    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: scl    (ADM_gall_1d_in,ADM_gall_1d ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: scl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!    real(RP), intent(in)  :: kh     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: kh     (ADM_gall_1d_in,ADM_gall_1d ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: kh_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: mfact

!    real(RP)  :: vxt    (ADM_gall   ,TI:TJ)
!    real(RP)  :: vyt    (ADM_gall   ,TI:TJ)
!    real(RP)  :: vzt    (ADM_gall   ,TI:TJ)
!    real(RP)  :: flux   (ADM_gall   ,AI:AJ)
    real(RP)  :: vxt    (ADM_gall_1d_in,ADM_gall_1d   ,TI:TJ+1)
    real(RP)  :: vyt    (ADM_gall_1d_in,ADM_gall_1d   ,TI:TJ+1)
    real(RP)  :: vzt    (ADM_gall_1d_in,ADM_gall_1d   ,TI:TJ+1)
    real(RP)  :: flux   (ADM_gall_1d_in,ADM_gall_1d   ,AI:AJ+1)

    real(RP) :: u1, u2, u3, smean

    integer :: nstart1, nend
    integer :: nstart2, nstart3
    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: n, k, l, v

    integer :: suf,i,j
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !fj add >
    OPRT_nstart = suf(ADM_gmin,ADM_gmin)
    OPRT_nend   = suf(ADM_gmax,ADM_gmax)
    !fj add<
    !---------------------------------------------------------------------------
    !cx  write(0,'(a,10i8)') "<OPRT_diffusion>"

!    call DEBUG_rapstart('OPRT_diffusion')
!    call start_collection('OPRT_diffusion')

       do l = 1, ADM_lall
!       nstart1 = suf(ADM_gmin-1,ADM_gmin-1)
!       nstart2 = suf(ADM_gmin-1,ADM_gmin)
!       nstart3 = suf(ADM_gmin,ADM_gmin-1)
!       nend   = suf(ADM_gmax  ,ADM_gmax  )

!$omp  parallel default(none),private(k,smean,u1,u2,&
!$omp  u3,i,j, vxt, vyt, vzt, flux),&
!$omp  shared(scl,ADM_gall_1d,ADM_gall_1d_in,ADM_gmax,&
!$omp  cinterp_TN,cinterp_TRA,kh,l,cinterp_HN,ADM_gmin, &
!$omp  ADM_have_sgp,cinterp_PRA,mfact,OPRT_nstart,dscl,OPRT_nend,ADM_kall)
!$omp do
       do k = 1, ADM_kall
          do j = ADM_gmin-1, ADM_gmax
#ifdef SINGLE_SIM
          do i = ADM_gmin-1, ADM_gmax/2
#else
          do i = ADM_gmin-1, ADM_gmax
#endif
!!OCL PREFETCH_READ(scl(i,j,k+1,l),level=2,strong=1)
          smean = ( scl(i,j,k,l) + scl(i+1,j,k,l) + scl(i+1,j+1,k,l) ) / 3._RP
          u1 = 0.5_RP * (scl(i,j    ,k,l)+scl(i+1,j  ,k,l)) - smean
          u2 = 0.5_RP * (scl(i+1,j  ,k,l)+scl(i+1,j+1,k,l)) - smean
          u3 = 0.5_RP * (scl(i+1,j+1,k,l)+scl(i,j    ,k,l)) - smean
          vxt(i,j,TI) = ( - u1 * cinterp_TN(i,j  ,l,AI ,1) &
                        - u2 * cinterp_TN(i+1,j,l,AJ ,1) &
                        + u3 * cinterp_TN(i,j  ,l,AIJ,1) ) * cinterp_TRA(i,j,l,TI)
          vyt(i,j,TI) = ( - u1 * cinterp_TN(i,j  ,l,AI ,2) &
                        - u2 * cinterp_TN(i+1,j,l,AJ ,2) &
                        + u3 * cinterp_TN(i,j  ,l,AIJ,2) ) * cinterp_TRA(i,j,l,TI)
          vzt(i,j,TI) = ( - u1 * cinterp_TN(i,j  ,l,AI ,3) &
                          - u2 * cinterp_TN(i+1,j,l,AJ ,3) &
                          + u3 * cinterp_TN(i,j  ,l,AIJ,3) ) * cinterp_TRA(i,j,l,TI)

          smean = ( scl(i,j,k,l) + scl(i+1,j+1,k,l) + scl(i,j+1,k,l) ) / 3.0_RP
          u1 = 0.5_RP * (scl(i,j    ,k,l)+scl(i+1,j+1,k,l)) - smean
          u2 = 0.5_RP * (scl(i+1,j+1,k,l)+scl(i,j+1  ,k,l)) - smean
          u3 = 0.5_RP * (scl(i,j+1  ,k,l)+scl(i,j    ,k,l)) - smean
          vxt(i,j,TJ) = ( - u1 * cinterp_TN(i,j  ,l,AIJ,1) &
                        + u2 * cinterp_TN(i,j+1,l,AI ,1) &
                        + u3 * cinterp_TN(i,j  ,l,AJ ,1) ) * cinterp_TRA(i,j,l,TJ)
          vyt(i,j,TJ) = ( - u1 * cinterp_TN(i,j  ,l,AIJ,2) &
                        + u2 * cinterp_TN(i,j+1,l,AI ,2) &
                        + u3 * cinterp_TN(i,j  ,l,AJ ,2) ) * cinterp_TRA(i,j,l,TJ)
          vzt(i,j,TJ) = ( - u1 * cinterp_TN(i,j  ,l,AIJ,3) &
                          + u2 * cinterp_TN(i,j+1,l,AI ,3) &
                          + u3 * cinterp_TN(i,j  ,l,AJ ,3) ) * cinterp_TRA(i,j,l,TJ)
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
!!$omp master
             vxt(ADM_gmin-1,ADM_gmin-1,TI) = vxt(ADM_gmin,ADM_gmin-1,TJ)
             vyt(ADM_gmin-1,ADM_gmin-1,TI) = vyt(ADM_gmin,ADM_gmin-1,TJ)
             vzt(ADM_gmin-1,ADM_gmin-1,TI) = vzt(ADM_gmin,ADM_gmin-1,TJ)
!!$omp end master
!!$omp barrier
       endif

       j = ADM_gmin-1
!!$omp do
#ifdef SINGLE_SIM
       do i = ADM_gmin-1, ADM_gmax/2
#else
       do i = ADM_gmin-1, ADM_gmax
#endif
          flux(i,j,AIJ) =   0.25_RP * ( (vxt(i,j  ,TI)+vxt(i,j  ,TJ)) * cinterp_HN(i,j,l,AIJ,1) &
                                   + (vyt(i,j  ,TI)+vyt(i,j  ,TJ)) * cinterp_HN(i,j,l,AIJ,2) &
                                   + (vzt(i,j  ,TI)+vzt(i,j  ,TJ)) * cinterp_HN(i,j,l,AIJ,3) &
                                   ) * (kh(i,j,k,l)+kh(i+1,j+1,k,l))
       enddo
!!$omp enddo nowait
!!$omp do
#ifdef SINGLE_SIM
       do i = ADM_gmin, ADM_gmax/2
#else
       do i = ADM_gmin, ADM_gmax
#endif
          flux(i,j,AJ ) =   0.25_RP * ( (vxt(i,j  ,TJ)+vxt(i-1,j,TI)) * cinterp_HN(i,j,l,AJ ,1) &
                                   + (vyt(i,j  ,TJ)+vyt(i-1,j,TI)) * cinterp_HN(i,j,l,AJ ,2) &
                                   + (vzt(i,j  ,TJ)+vzt(i-1,j,TI)) * cinterp_HN(i,j,l,AJ ,3) &
                                   ) * (kh(i,j,k,l)+kh(i,j+1  ,k,l))
       enddo
!!$omp enddo nowait

!!$omp do
       do j = ADM_gmin, ADM_gmax
#ifdef SINGLE_SIM
       do i = ADM_gmin-1, ADM_gmax /2
#else
       do i = ADM_gmin-1, ADM_gmax
#endif
!!OCL PREFETCH_READ(kh(i,j  ,k+1,l),level=2,strong=1)
!!OCL PREFETCH_READ(kh(i,j+1  ,k+1,l),level=2,strong=1)
          flux(i,j,AI ) = 0.25_RP *   ( (vxt(i,j-1,TJ)+vxt(i,j  ,TI)) * cinterp_HN(i,j,l,AI,1) &
                                   + (vyt(i,j-1,TJ)+vyt(i,j  ,TI)) * cinterp_HN(i,j,l,AI,2) &
                                   + (vzt(i,j-1,TJ)+vzt(i,j  ,TI)) * cinterp_HN(i,j,l,AI,3) &
                                   ) * (kh(i,j,k,l)+kh(i+1,j  ,k,l))
          flux(i,j,AIJ) =   0.25_RP * ( (vxt(i,j  ,TI)+vxt(i,j  ,TJ)) * cinterp_HN(i,j,l,AIJ,1) &
                                   + (vyt(i,j  ,TI)+vyt(i,j  ,TJ)) * cinterp_HN(i,j,l,AIJ,2) &
                                   + (vzt(i,j  ,TI)+vzt(i,j  ,TJ)) * cinterp_HN(i,j,l,AIJ,3) &
                                   ) * (kh(i,j,k,l)+kh(i+1,j+1,k,l))
       enddo
!       enddo
!!$omp enddo nowait
!!$omp do
!       do j = ADM_gmin-1, ADM_gmax
#ifdef SINGLE_SIM
       do i = ADM_gmin, ADM_gmax / 2 
#else
       do i = ADM_gmin, ADM_gmax
#endif
!!OCL PREFETCH_READ(kh(i,j  ,k+1,l),level=2,strong=1)
!!OCL PREFETCH_READ(kh(i,j+1  ,k+1,l),level=2,strong=1)
          flux(i,j,AJ ) =   0.25_RP * ( (vxt(i,j  ,TJ)+vxt(i-1,j,TI)) * cinterp_HN(i,j,l,AJ ,1) &
                                   + (vyt(i,j  ,TJ)+vyt(i-1,j,TI)) * cinterp_HN(i,j,l,AJ ,2) &
                                   + (vzt(i,j  ,TJ)+vzt(i-1,j,TI)) * cinterp_HN(i,j,l,AJ ,3) &
                                   ) * (kh(i,j,k,l)+kh(i,j+1  ,k,l))
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
!!$omp master
             flux(ADM_gmin,ADM_gmin-1,AJ) = 0._RP
!!$omp end master
!!$omp barrier
       endif

!!$omp do
       do j = ADM_gmin, ADM_gmax
#ifdef SINGLE_SIM
       do i = ADM_gmin, ADM_gmax/2
#else
       do i = ADM_gmin, ADM_gmax
#endif
!!OCL PREFETCH_WRITE(dscl(i,j,k+1,l),level=2)
          dscl(i,j,k,l) = ( flux(i,j,AI ) - flux(i-1,j  ,AI ) &
                        + flux(i,j,AIJ) - flux(i-1,j-1,AIJ) &
                        + flux(i,j,AJ ) - flux(i,j-1  ,AJ ) ) * cinterp_PRA(i,j,l) * mfact
       enddo
       enddo
!!$omp enddo nowait
       enddo
!$omp end parallel

    enddo

!cx cut off other loops

!    call DEBUG_rapend('OPRT_diffusion')
!    call stop_collection('OPRT_diffusion')

    return
  end subroutine OPRT_diffusion
end module mod_oprt

  subroutine dynamics_step
    use mod_precision
    use mod_adm, only: &
       ADM_gmax,    &
       ADM_gmin,    &
       ADM_gall_1d,    &
       ADM_gall_1d_in,    &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    use mod_oprt, only : &
        OPRT_snap_read, &
        OPRT_diffusion
    implicit none
    integer :: i
    !fj>
    integer :: k,l
    !fj<

    real(RP) :: dscl   (ADM_gall_1d_in,ADM_gall_1d   ,ADM_kall,ADM_lall   )
    real(RP) :: dscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: scl    (ADM_gall_1d_in,ADM_gall_1d   ,ADM_kall,ADM_lall   )
    real(RP) :: scl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: kh     (ADM_gall_1d_in,ADM_gall_1d   ,ADM_kall,ADM_lall   )
    real(RP) :: kh_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: mfact = 1.0_RP
    !fj>
    real(RP) :: sum1
    !fj<
#ifdef OUTPUTCHECK
    real(8) :: tmp_dscl = 0.0_8
    real(8) :: ref_dscl = 1.586542891295769E-03_8
    integer :: g,g1,g2
#endif

    PROF_INIT
    PROF_START_ALL

    call OPRT_snap_read( scl, scl_pl, kh, kh_pl, mfact )
!    write(0,'(a)') "starting the <OPRT_diffusion> loop"

!pre
!#ifdef SINGLE
!    call OPRT_diffusion( dscl, dscl_pl, scl, scl_pl, kh, kh_pl, mfact=1.0e0)
!#else
!    call OPRT_diffusion( dscl, dscl_pl, scl, scl_pl, kh, kh_pl, mfact=1.0d0)
!#endif
    !mfact= 1.0_RP
    !call OPRT_diffusion( dscl, dscl_pl, scl, scl_pl, kh, kh_pl, mfact)
    !cx call fapp_start( 'OPRT_diffusion', 1, 1 )
    PROF_START("OPRT_diffusion")
!    do i=1,60
    do i=1,1

    mfact= 1.0_RP
    call OPRT_diffusion( dscl, dscl_pl, scl, scl_pl, kh, kh_pl, mfact)
!#ifdef SINGLE
!    call OPRT_diffusion( dscl, dscl_pl, scl, scl_pl, kh, kh_pl, mfact=1.0e0)
!#else
!    call OPRT_diffusion( dscl, dscl_pl, scl, scl_pl, kh, kh_pl, mfact=1.0d0)
!#endif
    end do
    !cx call fapp_stop( 'OPRT_diffusion', 1, 1 )
    PROF_STOP("OPRT_diffusion")

    PROF_STOP_ALL
    PROF_FINALIZE

#ifdef OUTPUTCHECK
    do l=1, ADM_lall
     do k=1, ADM_kall
      do g1=ADM_gmin, ADM_gmax
#ifdef SINGLE_SIM
      do g2=ADM_gmin, ADM_gmax
#else
      do g2=ADM_gmin, ADM_gmax/2
#endif
        g = ADM_gall_1d * (g1 -1) + g2
          tmp_dscl = tmp_dscl + dble(dscl(g2,g1,k,l) ** 2)
      enddo
      enddo
     enddo
    enddo

!    write(*,*) "OUTPUTCHECK dscl=",sqrt(tmp_dscl)
    tmp_dscl = sqrt(tmp_dscl)
    call report_validation(tmp_dscl, ref_dscl, 0.00007_8)
#endif
!    write(0,'(a)') "finishing the <OPRT_diffusion> loop"
    return
  end subroutine dynamics_step

  subroutine get_walltime(wctime)
    use iso_fortran_env, only: int64
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp) :: wctime
    integer(int64) :: r, c
    call system_clock(c, r)
    wctime = real(c, dp) / r
  end subroutine get_walltime

program main
  ! use iso_fortran_env
  integer, parameter :: dp = kind(1.0d0)
  real(dp) :: start_time, end_time
  call get_walltime(start_time)
    call dynamics_step
  call get_walltime(end_time)
  write(0,'(a, F16.9)') "WALLTIME: ",end_time - start_time
end program main

