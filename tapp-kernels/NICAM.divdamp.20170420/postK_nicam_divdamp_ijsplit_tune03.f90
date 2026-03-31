#include "profiler.h"

program main
  integer, parameter :: dp = kind(1.0d0)
  real(dp) :: start_time, end_time
  call get_walltime(start_time)
    call dynamics_step
  call get_walltime(end_time)
  write(0,'(a, F16.9)') "WALLTIME: ",end_time - start_time
!   stop
end program main

!-------------------------------------------------------------------------------
!>
!! 3D Operator module
!!
!! @par Description
!!         This module contains the subroutines for differential oeprators using vertical metrics.
!!
!! @li  chosen for performance evaluation targetting post-K
!!
!<
module mod_adm
  implicit none

  integer, public  :: ADM_lall=1
  integer, public  :: ADM_lall_pl=2
  integer, public  :: ADM_gall=16900
  integer, public  :: ADM_gall_pl=6
!  integer, public  :: ADM_kall=96
  integer, public  :: ADM_kall=48
!  integer, public  :: ADM_gall_1d=130
  integer, public,parameter  :: ADM_gall_1d=130
#ifdef SINGLE_SIM
  integer, public,parameter  :: ADM_gall_1d_in=65
#else
  integer, public,parameter  :: ADM_gall_1d_in=130
#endif
  integer, public  :: ADM_gmin=2
  integer, public  :: ADM_gmax=129
  integer, public  :: ADM_gmax_pl=6
  integer, public  :: ADM_kmin=2
!  integer, public  :: ADM_kmax=95
  integer, public  :: ADM_kmax=47
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
module mod_grd
    use mod_precision
    use mod_adm
  real(RP), public, save ::  GRD_rdgz (ADM_kall_para)
end module mod_grd
module mod_vmtr
    use mod_precision
    use mod_adm
    real(RP), public :: VMTR_RGAM       (ADM_gall_para,ADM_kall_para,ADM_lall_para ) 
    real(RP), public :: VMTR_RGAMH      (ADM_gall_para,ADM_kall_para,ADM_lall_para ) 
    real(RP), public :: VMTR_RGSQRTH    (ADM_gall_para,ADM_kall_para,ADM_lall_para ) 
    real(RP), public :: VMTR_C2WfactGz(ADM_gall_para,6,ADM_kall_para,ADM_lall_para )

!    real(RP), public, target :: RGAM (ADM_gall_para,ADM_kall_para,ADM_lall_para )
!    real(RP), public, target :: RGAMH (ADM_gall_para,ADM_kall_para,ADM_lall_para )
!    real(RP), public, target :: RGSQRTH (ADM_gall_para,ADM_kall_para,ADM_lall_para )
!    real(RP), public, target :: C2WfactGz(6,ADM_gall_para,ADM_kall_para,ADM_lall_para )
!    real(RP), public, pointer :: VMTR_RGAM (:,:,:)
!    real(RP), public, pointer :: VMTR_RGAMH (:,:,:)
!    real(RP), public, pointer :: VMTR_RGSQRTH (:,:,:)
!    real(RP), public, pointer :: VMTR_C2WfactGz(:,:,:,:)

!    real(RP), public, target :: RGAM (ADM_gall_1d_para,ADM_gall_1d_para,ADM_kall_para,ADM_lall_para )
!    real(RP), public, target :: RGAMH (ADM_gall_1d_para,ADM_gall_1d_para,ADM_kall_para,ADM_lall_para )
!    real(RP), public, target :: RGSQRTH (ADM_gall_1d_para,ADM_gall_1d_para,ADM_kall_para,ADM_lall_para )
!    real(RP), public, target :: C2WfactGz(6,ADM_gall_1d_para,ADM_gall_1d_para,ADM_kall_para,ADM_lall_para )
!    real(RP), public, pointer :: VMTR_RGAM (:,:,:,:)
!    real(RP), public, pointer :: VMTR_RGAMH (:,:,:,:)
!    real(RP), public, pointer :: VMTR_RGSQRTH (:,:,:,:)
!    real(RP), public, pointer :: VMTR_C2WfactGz(:,:,:,:,:)

!    real(RP), public :: VMTR_RGAM       (ADM_gall_1d_para,ADM_gall_1d_para,ADM_kall_para,ADM_lall_para)
!    real(RP), public :: VMTR_RGAMH      (ADM_gall_1d_para,ADM_gall_1d_para,ADM_kall_para,ADM_lall_para )
!    real(RP), public :: VMTR_RGSQRTH    (ADM_gall_1d_para,ADM_gall_1d_para,ADM_kall_para,ADM_lall_para )
!    real(RP), public :: VMTR_C2WfactGz(6,ADM_gall_1d_para,ADM_gall_1d_para,ADM_kall_para,ADM_lall_para )

end module mod_vmtr
module mod_oprt
    use mod_precision
    use mod_adm
!  real(RP), public, target, save :: TN(ADM_AI:ADM_AJ,3,ADM_gall_para,ADM_lall_para)
!  real(RP), public, target, save :: HN(ADM_AI:ADM_AJ,3,ADM_gall_para,ADM_lall_para)
!  real(RP), public, target, save :: TRA(ADM_TI:ADM_TJ,ADM_gall_para,ADM_lall_para)
!  real(RP), public, target, save :: PRA(ADM_gall_para,ADM_lall_para)
!  real(RP), public, pointer, save ::  cinterp_TN(:,:,:,:)
!  real(RP), public, pointer, save ::  cinterp_HN(:,:,:,:)
!  real(RP), public, pointer, save :: cinterp_TRA(:,:,:)
!  real(RP), public, pointer, save :: cinterp_PRA(:,:)

  !real(RP), public, target, save :: TN(ADM_AI:ADM_AJ,3,ADM_gall_1d_para,ADM_gall_1d_para,ADM_lall_para)
  !real(RP), public, target, save :: HN(ADM_AI:ADM_AJ,3,ADM_gall_1d_para,ADM_gall_1d_para,ADM_lall_para)
!  real(RP), public, target, save :: TN(3,ADM_gall_1d_para,ADM_gall_1d_para,ADM_AI:ADM_AJ,ADM_lall_para)
!  real(RP), public, target, save :: HN(3,ADM_gall_1d_para,ADM_gall_1d_para,ADM_AI:ADM_AJ,ADM_lall_para)
!  real(RP), public, target, save :: TRA(ADM_TI:ADM_TJ,ADM_gall_1d_para,ADM_gall_1d_para,ADM_lall_para)
!  real(RP), public, target, save :: PRA(ADM_gall_1d_para,ADM_gall_1d_para,ADM_lall_para)
!  real(RP), public, save ::  cinterp_TN(ADM_gall_1d_para,ADM_gall_1d_para,3,ADM_AI:ADM_AJ,ADM_lall_para)
!  real(RP), public, save :: cinterp_TN(ADM_AI:ADM_AJ,3,ADM_gall_1d_para,ADM_gall_1d_para,ADM_lall_para)
!!  real(RP), public, save ::  cinterp_HN(ADM_gall_1d_para,ADM_gall_1d_para,3,ADM_AI:ADM_AJ,ADM_lall_para)
!  real(RP), public, save :: cinterp_HN(ADM_AI:ADM_AJ,3,ADM_gall_1d_para,ADM_gall_1d_para,ADM_lall_para)
!  real(RP), public, save :: cinterp_TRA(ADM_TI:ADM_TJ,ADM_gall_1d_para,ADM_gall_1d_para,ADM_lall_para)
!  real(RP), public, save :: cinterp_PRA(ADM_gall_1d_para,ADM_gall_1d_para,ADM_lall_para)
  real(RP), public, save ::  cinterp_TN(ADM_gall_para,ADM_AI:ADM_AJ,3,ADM_lall_para)
  real(RP), public, save ::  cinterp_HN(ADM_gall_para,ADM_AI:ADM_AJ,3,ADM_lall_para)
  real(RP), public, save :: cinterp_TRA(ADM_gall_para,ADM_TI:ADM_TJ,ADM_lall_para)
  real(RP), public, save :: cinterp_PRA(ADM_gall_para,ADM_lall_para)
end module mod_oprt
module mod_oprt3d
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !

    use mod_precision
    use mod_adm, only: &
       ADM_have_sgp, &
       ADM_lall,     &
       ADM_lall_pl,  &
       ADM_gall,     &
       ADM_gall_pl,  &
       ADM_kall,     &
       ADM_gall_1d,  &
       ADM_gall_1d_in,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_gmax_pl,  &
       ADM_kmin,     &
      ADM_kall_para, &
      ADM_gall_para, &
      ADM_lall_para, &
     TI  => ADM_TI,  &
     TJ  => ADM_TJ,  &
     AI  => ADM_AI,  &
     AIJ => ADM_AIJ, &
     AJ  => ADM_AJ,  &
       ADM_kmax
    use mod_grd, only: &
       GRD_rdgz
    use mod_vmtr, only: &
       VMTR_RGAM,        &
       VMTR_RGAMH,       &
       VMTR_RGSQRTH,     &
       VMTR_C2WfactGz
    use mod_oprt

  implicit none
  private

  !++ Public procedure
  public :: OPRT3D_snap_read
  public :: OPRT3D_divdamp
!  real(RP), public, save ::  cinterp_TN(AI:AJ,3,ADM_gall_para,ADM_lall_para)
!  real(RP), public, save ::  cinterp_HN(AI:AJ,3,ADM_gall_para,ADM_lall_para)
!  real(RP), public, save :: cinterp_TRA(TI:TJ,ADM_gall_para,ADM_lall_para)
!  real(RP), public, save :: cinterp_PRA(ADM_gall_para,ADM_lall_para)

  !++ Public parameters & variables
  !++ Private procedures
  !++ Private parameters & variables

contains

  subroutine OPRT3D_divdamp( ddivdx, ddivdy, ddivdz, &
                            rhogvx, rhogvy, rhogvz, rhogw,&
                            cinterp_TN_l,cinterp_HN_l,cinterp_TRA_l,cinterp_PRA_l, &
                            VMTR_RGAM_l,VMTR_RGAMH_l,VMTR_RGSQRTH_l,VMTR_C2WfactGz_l)



    implicit none
    real(RP), intent(out) :: ddivdx   (ADM_gall   ,ADM_kall,ADM_lall   ) ! tendency
    real(RP), intent(out) :: ddivdy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: ddivdz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvx   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*vx { gam2 x G^1/2 }
    real(RP), intent(in)  :: rhogvy   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*vy { gam2 x G^1/2 }
    real(RP), intent(in)  :: rhogvz   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*vz { gam2 x G^1/2 }
    real(RP), intent(in)  :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w  { gam2 x G^1/2 }

    real(RP) :: cinterp_TN_l(ADM_gall_para,ADM_AI:ADM_AJ,3,ADM_lall_para)
    real(RP) :: cinterp_HN_l(ADM_gall_para,ADM_AI:ADM_AJ,3,ADM_lall_para)
    real(RP) :: cinterp_TRA_l(ADM_gall_para,ADM_TI:ADM_TJ,ADM_lall_para) 
    real(RP) :: cinterp_PRA_l(ADM_gall_para,ADM_lall_para)
    real(RP) :: VMTR_RGAM_l(ADM_gall_para,ADM_kall_para,ADM_lall_para )
    real(RP) :: VMTR_RGAMH_l(ADM_gall_para,ADM_kall_para,ADM_lall_para )
    real(RP) :: VMTR_RGSQRTH_l(ADM_gall_para,ADM_kall_para,ADM_lall_para )
    real(RP) :: VMTR_C2WfactGz_l(ADM_gall_para,6,ADM_kall_para,ADM_lall_para )

!fj    real(RP) :: sclt         (ADM_gall   ,ADM_kall,TI:TJ) ! scalar on the hexagon vertex
    real(RP) :: sclt         (ADM_gall   ,TI:TJ+1) ! scalar on the hexagon vertex
    real(RP) :: sclt_rhogw

!fj    real(RP) :: rhogvx_vm   (ADM_gall   ,ADM_kall) ! rho*vx / vertical metrics
!fj    real(RP) :: rhogvy_vm   (ADM_gall   ,ADM_kall) ! rho*vy / vertical metrics
!fj    real(RP) :: rhogvz_vm   (ADM_gall   ,ADM_kall) ! rho*vz / vertical metrics
    real(RP) :: rhogvx_vm   (ADM_gall) ! rho*vx / vertical metrics
    real(RP) :: rhogvy_vm   (ADM_gall) ! rho*vy / vertical metrics
    real(RP) :: rhogvz_vm   (ADM_gall) ! rho*vz / vertical metrics
    real(RP) :: rhogw_vm    (ADM_gall,   ADM_kall) ! rho*w  / vertical metrics

    integer :: nstart, nend
    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: g, k, l, v, n

    integer :: suf,i,j
    integer :: OPRT_nstart,OPRT_nend
    real(RP) :: x1,y1,z1,x2,y2,z2
    real(RP) :: wk1,wk2,wk3,wk4,wk5,wk6
    suf(i,j) = ADM_gall_1d * ((j)-1) + (i)
    !---------------------------------------------------------------------------
    OPRT_nstart = suf(ADM_gmin,ADM_gmin)
    OPRT_nend   = suf(ADM_gmax,ADM_gmax)
   


    
    do l = 1, ADM_lall
!$omp parallel do private(k,g)
       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          rhogw_vm(g,k) = ( VMTR_C2WfactGz_l(g,1,k,l) * rhogvx(g,k  ,l) &
                          + VMTR_C2WfactGz_l(g,2,k,l) * rhogvx(g,k-1,l) &
                          + VMTR_C2WfactGz_l(g,3,k,l) * rhogvy(g,k  ,l) &
                          + VMTR_C2WfactGz_l(g,4,k,l) * rhogvy(g,k-1,l) &
                          + VMTR_C2WfactGz_l(g,5,k,l) * rhogvz(g,k  ,l) &
                          + VMTR_C2WfactGz_l(g,6,k,l) * rhogvz(g,k-1,l) &
                          ) * VMTR_RGAMH_l(g,k,l)                       & ! horizontal contribution
                        + rhogw(g,k,l) * VMTR_RGSQRTH_l(g,k,l)            ! vertical   contribution
       enddo

       enddo
!$omp parallel do private(g)
       do g = 1, ADM_gall
          rhogw_vm(g,ADM_kmin  ) = 0.0_RP
          rhogw_vm(g,ADM_kmax+1) = 0.0_RP
       enddo
!fj>
       nstart = suf(ADM_gmin-1,ADM_gmin-1)
       nend   = suf(ADM_gmax,  ADM_gmax  )
!fj<
!$omp parallel do private(rhogvx_vm,rhogvy_vm,rhogvz_vm,ij,ip1j,ip1jp1,sclt_rhogw,ijp1,im1j,sclt,ijm1,im1jm1, &
!$omp                     x1,y1,z1,x2,y2,z2,wk1,wk2,wk3,wk4,wk5,wk6)
       do k = ADM_kmin, ADM_kmax
       do g = 1, ADM_gall
!fj          rhogvx_vm(g,k) = rhogvx(g,k,l) * VMTR_RGAM(g,k,l)
!fj          rhogvy_vm(g,k) = rhogvy(g,k,l) * VMTR_RGAM(g,k,l)
!fj          rhogvz_vm(g,k) = rhogvz(g,k,l) * VMTR_RGAM(g,k,l)
          rhogvx_vm(g) = rhogvx(g,k,l) * VMTR_RGAM_l(g,k,l)
          rhogvy_vm(g) = rhogvy(g,k,l) * VMTR_RGAM_l(g,k,l)
          rhogvz_vm(g) = rhogvz(g,k,l) * VMTR_RGAM_l(g,k,l)
       enddo
!fj       enddo

!fj       nstart = suf(ADM_gmin-1,ADM_gmin-1)
!fj       nend   = suf(ADM_gmax,  ADM_gmax  )

!fj       do k = ADM_kmin, ADM_kmax
!ocl loop_nofusion
       do ij = nstart, nend
          x1 = - (rhogvx_vm(ij              )+rhogvx_vm(ij+1            )) * cinterp_TN_l(ij,               AI ,1,l) &
               - (rhogvx_vm(ij+1            )+rhogvx_vm(ij+1+ADM_gall_1d)) * cinterp_TN_l(ij+1,             AJ ,1,l) &
               + (rhogvx_vm(ij+1+ADM_gall_1d)+rhogvx_vm(ij              )) * cinterp_TN_l(ij,               AIJ,1,l) 
          x2 = - (rhogvx_vm(ij              )+rhogvx_vm(ij+1+ADM_gall_1d)) * cinterp_TN_l(ij,               AIJ,1,l) &
               + (rhogvx_vm(ij+1+ADM_gall_1d)+rhogvx_vm(ij  +ADM_gall_1d)) * cinterp_TN_l(ij   +ADM_gall_1d,AI ,1,l) &
               + (rhogvx_vm(ij  +ADM_gall_1d)+rhogvx_vm(ij              )) * cinterp_TN_l(ij,               AJ ,1,l) 

          y1 = - (rhogvy_vm(ij              )+rhogvy_vm(ij+1            )) * cinterp_TN_l(ij,               AI ,2,l) &
               - (rhogvy_vm(ij+1            )+rhogvy_vm(ij+1+ADM_gall_1d)) * cinterp_TN_l(ij+1,             AJ ,2,l) &
               + (rhogvy_vm(ij+1+ADM_gall_1d)+rhogvy_vm(ij              )) * cinterp_TN_l(ij,               AIJ,2,l) 
          y2 = - (rhogvy_vm(ij              )+rhogvy_vm(ij+1+ADM_gall_1d)) * cinterp_TN_l(ij,               AIJ,2,l) &
               + (rhogvy_vm(ij+1+ADM_gall_1d)+rhogvy_vm(ij  +ADM_gall_1d)) * cinterp_TN_l(ij   +ADM_gall_1d,AI ,2,l) &
               + (rhogvy_vm(ij  +ADM_gall_1d)+rhogvy_vm(ij              )) * cinterp_TN_l(ij,               AJ ,2,l) 

          z1 = - (rhogvz_vm(ij              )+rhogvz_vm(ij+1            )) * cinterp_TN_l(ij,               AI ,3,l) &
               - (rhogvz_vm(ij+1            )+rhogvz_vm(ij+1+ADM_gall_1d)) * cinterp_TN_l(ij+1,             AJ ,3,l) &
               + (rhogvz_vm(ij+1+ADM_gall_1d)+rhogvz_vm(ij              )) * cinterp_TN_l(ij,               AIJ,3,l) 
          z2 = - (rhogvz_vm(ij              )+rhogvz_vm(ij+1+ADM_gall_1d)) * cinterp_TN_l(ij,               AIJ,3,l) &
               + (rhogvz_vm(ij+1+ADM_gall_1d)+rhogvz_vm(ij  +ADM_gall_1d)) * cinterp_TN_l(ij   +ADM_gall_1d,AI ,3,l) &
               + (rhogvz_vm(ij  +ADM_gall_1d)+rhogvz_vm(ij              )) * cinterp_TN_l(ij,               AJ ,3,l) 

          sclt(ij,TI) = ( x1 + y1 + z1 ) * 0.5_RP * cinterp_TRA_l(ij,TI,l)
          sclt(ij,TJ) = ( x2 + y2 + z2 ) * 0.5_RP * cinterp_TRA_l(ij,TJ,l) 
       enddo

!ocl loop_nofusion
       do ij = nstart, nend
          sclt(ij,TI) = ( ( rhogw_vm(ij,k+1) + rhogw_vm(ij+1,k+1) + rhogw_vm(ij+1+ADM_gall_1d,k+1) ) &
                        - ( rhogw_vm(ij,k  ) + rhogw_vm(ij+1,k  ) + rhogw_vm(ij+1+ADM_gall_1d,k  ) ) &
                        ) / 3.0_RP * GRD_rdgz(k) + sclt(ij,TI)
          sclt(ij,TJ) = ( ( rhogw_vm(ij,k+1) + rhogw_vm(ij+ADM_gall_1d,k+1) + rhogw_vm(ij+1+ADM_gall_1d,k+1) ) &
                        - ( rhogw_vm(ij,k  ) + rhogw_vm(ij+ADM_gall_1d,k  ) + rhogw_vm(ij+1+ADM_gall_1d,k  ) ) &
                        ) / 3.0_RP * GRD_rdgz(k) + sclt(ij,TJ)
       enddo

       do ij = OPRT_nstart, OPRT_nend

          wk1 = sclt(ij  -ADM_gall_1d,TJ) + sclt(ij,              TI)
          wk2 = sclt(ij,              TI) + sclt(ij,              TJ)
          wk3 = sclt(ij,              TJ) + sclt(ij-1,            TI)
          wk4 = sclt(ij-1-ADM_gall_1d,TJ) + sclt(ij-1,            TI)
          wk5 = sclt(ij-1-ADM_gall_1d,TI) + sclt(ij-1-ADM_gall_1d,TJ)
          wk6 = sclt(ij  -ADM_gall_1d,TJ) + sclt(ij-1-ADM_gall_1d,TI)

          ddivdx(ij,k,l) = (  wk1 * cinterp_HN_l(ij,    AI ,1,l) &
                            + wk2 * cinterp_HN_l(ij,    AIJ,1,l) &
                            + wk3 * cinterp_HN_l(ij,    AJ ,1,l) &
                            - wk4 * cinterp_HN_l(ij-1,  AI ,1,l) &
                            - wk5 * cinterp_HN_l(ij-1-ADM_gall_1d,AIJ,1,l) &
                            - wk6 * cinterp_HN_l(ij  -ADM_gall_1d,AJ ,1,l) &
                          ) !* 0.5_RP * cinterp_PRA_l(ij,l)                               
                                                                                      
          ddivdy(ij,k,l) = (  wk1 * cinterp_HN_l(ij,    AI ,2,l) &
                            + wk2 * cinterp_HN_l(ij,    AIJ,2,l) &
                            + wk3 * cinterp_HN_l(ij,    AJ ,2,l) &
                            - wk4 * cinterp_HN_l(ij-1,  AI ,2,l) &
                            - wk5 * cinterp_HN_l(ij-1-ADM_gall_1d,AIJ,2,l) &
                            - wk6 * cinterp_HN_l(ij  -ADM_gall_1d,AJ ,2,l) &
                          ) !* 0.5_RP * cinterp_PRA_l(ij,l)

          ddivdz(ij,k,l) = (  wk1 * cinterp_HN_l(ij,    AI ,3,l) &
                            + wk2 * cinterp_HN_l(ij,    AIJ,3,l) &
                            + wk3 * cinterp_HN_l(ij,    AJ ,3,l) &
                            - wk4 * cinterp_HN_l(ij-1,  AI ,3,l) &
                            - wk5 * cinterp_HN_l(ij-1-ADM_gall_1d,AIJ,3,l) &
                            - wk6 * cinterp_HN_l(ij  -ADM_gall_1d,  AJ ,3,l) &
                          ) !* 0.5_RP * cinterp_PRA_l(ij,l)
       end do
!ocl loop_nofusion
       do ij = OPRT_nstart, OPRT_nend
          ddivdx(ij,k,l) = ddivdx(ij,k,l) * 0.5_RP * cinterp_PRA_l(ij,l)
          ddivdy(ij,k,l) = ddivdy(ij,k,l) * 0.5_RP * cinterp_PRA_l(ij,l)
          ddivdz(ij,k,l) = ddivdz(ij,k,l) * 0.5_RP * cinterp_PRA_l(ij,l)
       enddo
       enddo

!$omp parallel do private(g)
       do g = 1, ADM_gall
          ddivdx(g,ADM_kmin-1,l) = 0.0_RP
          ddivdx(g,ADM_kmax+1,l) = 0.0_RP
          ddivdy(g,ADM_kmin-1,l) = 0.0_RP
          ddivdy(g,ADM_kmax+1,l) = 0.0_RP
          ddivdz(g,ADM_kmin-1,l) = 0.0_RP
          ddivdz(g,ADM_kmax+1,l) = 0.0_RP
       enddo
    enddo

    return
  end subroutine OPRT3D_divdamp


  subroutine OPRT3D_snap_read( rhogvx, rhogvy, rhogvz, rhogw )

    implicit none
    real(RP), intent(out)  :: rhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out)  :: rhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out)  :: rhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out)  :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
!    real(RP), intent(out)  :: rhogvx   (ADM_gall_1d_in,ADM_gall_1d   ,ADM_kall,ADM_lall   )
!    real(RP), intent(out)  :: rhogvy   (ADM_gall_1d_in,ADM_gall_1d   ,ADM_kall,ADM_lall   )
!    real(RP), intent(out)  :: rhogvz   (ADM_gall_1d_in,ADM_gall_1d   ,ADM_kall,ADM_lall   )
!    real(RP), intent(out)  :: rhogw    (ADM_gall_1d_in,ADM_gall_1d   ,ADM_kall,ADM_lall   )
    integer :: i,j,k,l,t,g

!cx Better To Do: first touch the following arrays


#ifndef OUTPUTCHECK
        cinterp_TN( 1:ADM_gall,AI:AJ,1:3,1:ADM_lall) = 1.0
        cinterp_HN( 1:ADM_gall,AI:AJ,1:3,1:ADM_lall) = 1.0
        cinterp_TRA(1:ADM_gall,TI:TJ,    1:ADM_lall) = 1.0
        cinterp_PRA(1:ADM_gall,          1:ADM_lall) = 1.0
        rhogvx (1:ADM_gall ,1:ADM_kall, 1:ADM_lall) = 1.0
        rhogvy (1:ADM_gall ,1:ADM_kall, 1:ADM_lall) = 1.0
        rhogvz (1:ADM_gall ,1:ADM_kall, 1:ADM_lall) = 1.0
        rhogw  (1:ADM_gall ,1:ADM_kall, 1:ADM_lall) = 1.0

!        cinterp_TN(AI:AJ,1:3,1:ADM_gall_1d,1:ADM_gall_1d,1:ADM_lall) = 1.0
!        cinterp_HN(AI:AJ,1:3,1:ADM_gall_1d,1:ADM_gall_1d,1:ADM_lall) = 1.0
!        cinterp_TN(1:ADM_gall_1d_in,1:ADM_gall_1d,1:3,AI:AJ,1:ADM_lall) = 1.0
!        cinterp_HN(1:ADM_gall_1d_in,1:ADM_gall_1d,1:3,AI:AJ,1:ADM_lall) = 1.0
!        cinterp_TRA(TI:TJ,1:ADM_gall_1d,1:ADM_gall_1d,1:ADM_lall) = 1.0
!        cinterp_TRA(TI:TJ,1:ADM_gall_1d_in,1:ADM_gall_1d,1:ADM_lall) = 1.0
!        cinterp_PRA(      1:ADM_gall_1d_in,1:ADM_gall_1d,1:ADM_lall) = 1.0
!        rhogvx (1:ADM_gall_1d_in,1:ADM_gall_1d ,1:ADM_kall, 1:ADM_lall) = 1.0
!        rhogvy (1:ADM_gall_1d_in,1:ADM_gall_1d ,1:ADM_kall, 1:ADM_lall) = 1.0
!        rhogvz (1:ADM_gall_1d_in,1:ADM_gall_1d ,1:ADM_kall, 1:ADM_lall) = 1.0
!        rhogw  (1:ADM_gall_1d_in,1:ADM_gall_1d ,1:ADM_kall, 1:ADM_lall) = 1.0

        VMTR_C2WfactGz(:,:,:,:)=1.0
        VMTR_RGAM(:,:,:)=1.0
        VMTR_RGAMH(:,:,:)=1.0
        VMTR_RGSQRTH(:,:,:)=1.0

!        VMTR_C2WfactGz(:,:,:,:,:)=1.0
!        VMTR_RGAM(:,:,:,:)=1.0
!        VMTR_RGAMH(:,:,:,:)=1.0
!        VMTR_RGSQRTH(:,:,:,:)=1.0

        GRD_rdgz(:)=1.0
#else

!ocl serial
        do t = 1,3
         do i = 1, ADM_gall_1d
          do j = 1, ADM_gall_1d_in
          g = ADM_gall_1d * (j -1) + i
           do k = AI, AJ
            do l = 1, ADM_lall
              cinterp_TN(g,k,t,l) = 0.1 / t * i * j / k * l / 1000
              cinterp_HN(g,k,t,l) = 0.2 / t * i * j / k * l / 1000
            enddo
           enddo
          enddo
         enddo
        enddo
!ocl serial
        do t = TI,TJ
         do i = 1, ADM_gall_1d
          do j = 1, ADM_gall_1d_in
          g = ADM_gall_1d * (j -1) + i
           do l = 1, ADM_lall
              cinterp_TRA(g,t,l) = 0.03 / t * i * j * l / 1000
           enddo
          enddo
         enddo
        enddo
!ocl serial
        do i = 1, ADM_gall_1d
         do j = 1, ADM_gall_1d_in
          g = ADM_gall_1d * (j -1) + i
          do l = 1, ADM_lall
            cinterp_PRA(g,l) = 0.04 * i * j *  l / 1000
          enddo
         enddo
        enddo

!ocl serial
        do i = 1, ADM_gall_1d
         do j = 1, ADM_gall_1d_in
          g = ADM_gall_1d * (j -1) + i
          do k = 1, ADM_kall
           do l = 1, ADM_lall
             rhogvx (g,k,l) = 0.05 * i * j / k * l / 1000
             rhogvy (g,k,l) = 0.06 * i * j / k * l / 1000
             rhogvz (g,k,l) = 0.07 * i * j / k * l / 1000
             rhogw  (g,k,l) = 0.08 * i * j / k * l / 1000
           enddo
          enddo
         enddo
        enddo

!ocl serial
        do t = 1,6
         do i = 1, ADM_gall_1d
          do j = 1, ADM_gall_1d_in
           g = ADM_gall_1d * (j -1) + i
           do k = 1, ADM_kall
            do l = 1, ADM_lall
            !!!  VMTR_C2WfactGz(t,g,k,l) = 0.0001 * i*j / k*l / t / 1000
              VMTR_C2WfactGz(g,t,k,l) = 0.0001 * i*j / k*l / t / 1000
            enddo
           enddo
          enddo
         enddo
        enddo

!ocl serial
        do i = 1, ADM_gall_1d
         do j = 1, ADM_gall_1d_in
          g = ADM_gall_1d * (j -1) + i
          do k = 1, ADM_kall
           do l = 1, ADM_lall
             VMTR_RGAM(g,k,l)= 0.0003 * i*j/k*l / 1000
             VMTR_RGAMH(g,k,l)= 0.0004 * i*j/k*l / 1000
             VMTR_RGSQRTH(g,k,l)=  0.0005 * i*j/k*l / 1000
           enddo
          enddo
         enddo
        enddo

!ocl serial
        do k = 1, ADM_kall
         GRD_rdgz(k)=k * 1.0
        enddo

#endif
        return
  end subroutine OPRT3D_snap_read



end module mod_oprt3d
!-------------------------------------------------------------------------------


  subroutine dynamics_step
    use mod_precision
    use mod_oprt3d, only : &
        OPRT3D_snap_read, &
        OPRT3D_divdamp
    use mod_adm, only: &
       ADM_gall_1d,    &
       ADM_gall_1d_in,    &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_kall
    use mod_oprt
    use mod_vmtr, only: &
       VMTR_RGAM,        &
       VMTR_RGAMH,       &
       VMTR_RGSQRTH,     &
       VMTR_C2WfactGz
    implicit none
    integer :: i
    !fj<

    real(RP) :: ddivdx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ddivdy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: ddivdz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: ddivdz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
!cx Better To Do: first touch above arrays
#ifdef USE_TIMER
    real*8 t_all_0, t_all, t_kernel_0, t_kernel
#endif
#if OUTPUTCHECK
    real(8) :: tmp_ddivdx, tmp_ddivdy, tmp_ddivdz
    integer :: l,k,g,g1,g2
!    real(8) :: ref_ddivdx = 6.142981652161969E-03_8
!    real(8) :: ref_ddivdy = 3.071490826080985E-03_8
!    real(8) :: ref_ddivdz = 2.047660499908134E-03_8
!    real(8) :: ref_ddivdx = 6.155367326516758E-03_8
!    real(8) :: ref_ddivdy = 3.077683663258379E-03_8
!    real(8) :: ref_ddivdz = 2.051789127012335E-03_8
    real(8) :: ref_ddivdx = 6.353791446821666E-03_8
    real(8) :: ref_ddivdy = 3.176895723410833E-03_8
    real(8) :: ref_ddivdz = 2.117930394569597E-03_8
#endif
    PROF_INIT
    PROF_START_ALL
#ifdef USE_TIMER
    call gettod(t_all_0)
#endif

!    call org_allocate()

    call OPRT3D_snap_read( rhogvx, rhogvy, rhogvz, rhogw )
#ifndef SIM_MEASURE
    call OPRT3D_divdamp( ddivdx, ddivdy, ddivdz, &
                            rhogvx, rhogvy, rhogvz, rhogw  , cinterp_TN, cinterp_HN,cinterp_TRA,cinterp_PRA, &
                            VMTR_RGAM,VMTR_RGAMH,VMTR_RGSQRTH,VMTR_C2WfactGz)
#endif
    PROF_START("divdamp")
#ifdef USE_TIMER
    call gettod(t_kernel_0)
#endif

!    do i=1,140
    do i=1,1
    call OPRT3D_divdamp( ddivdx, ddivdy, ddivdz, &
                            rhogvx, rhogvy, rhogvz, rhogw  , cinterp_TN, cinterp_HN,cinterp_TRA,cinterp_PRA, &
                            VMTR_RGAM,VMTR_RGAMH,VMTR_RGSQRTH,VMTR_C2WfactGz)
    end do
    PROF_STOP("divdamp")
#ifdef USE_TIMER
    call gettod(t_kernel)
    t_kernel = t_kernel - t_kernel_0
#endif
#ifdef USE_TIMER
    call gettod(t_all)
    t_all = t_all - t_all_0
#endif
    PROF_STOP_ALL
    PROF_FINALIZE
#ifdef USE_TIMER
    write(*,*)"t_kernel ",t_kernel
    write(*,*)"t_all ",t_all
#endif

#ifdef OUTPUTCHECK
    
    tmp_ddivdx = 0.0
    tmp_ddivdy = 0.0
    tmp_ddivdz = 0.0
!ocl serial
    do l=1, ADM_lall
     do k=1, ADM_kall
!      do g=1, ADM_gall
      do g1=ADM_gmin, ADM_gmax
#ifdef SINGLE_SIM
      do g2=ADM_gmin, ADM_gmax/2
#else
      do g2=ADM_gmin, ADM_gmax
#endif
        g = ADM_gall_1d * (g1 -1) + g2
        tmp_ddivdx = tmp_ddivdx + dble((ddivdx(g,k,l) ** 2))
        tmp_ddivdy = tmp_ddivdy + dble((ddivdy(g,k,l) ** 2))
        tmp_ddivdz = tmp_ddivdz + dble((ddivdz(g,k,l) ** 2))
      enddo
      enddo
     enddo
    enddo

 
!    call report_validation(sqrt(tmp_ddivdx),ref_ddivdx,0.1)
!    call report_validation(sqrt(tmp_ddivdy),ref_ddivdy,0.1)
!    call report_validation(sqrt(tmp_ddivdz),ref_ddivdz,0.1)

    tmp_ddivdx = sqrt(tmp_ddivdx)
    tmp_ddivdy = sqrt(tmp_ddivdy)
    tmp_ddivdz = sqrt(tmp_ddivdz)

    call report_validation(tmp_ddivdx,ref_ddivdx,0.00003_8)
    call report_validation(tmp_ddivdy,ref_ddivdy,0.00003_8)
    call report_validation(tmp_ddivdz,ref_ddivdz,0.00003_8)

#endif

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

