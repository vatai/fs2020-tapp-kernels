#include "profiler.h"

program main
    call dynamics_step

!FJ<
!   stop
end program main


  subroutine dynamics_step
    use mod_precision
    use mod_vi, only : &
        vi_rhow_solver
    use mod_adm, only: &
        ADM_gall,    &
        ADM_gall_1d,    &
        ADM_gall_1d_in,    &
        ADM_gall_pl, &
        ADM_lall,    &
        ADM_lall_pl, &
        ADM_kall
    use mod_vmtr
    use mod_vi
    use mod_grd 

    implicit none
    integer :: i
    real(RP) rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) rhogw0_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) preg0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) rhog0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) Sr_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) Sw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) Sp_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) dt

    integer g,k,l,j
    real(RP) rhogw_ij    (ADM_gall_1d_in,ADM_kall,ADM_gall_1d,ADM_lall   )
    real(RP) rhogw0_ij   (ADM_gall_1d_in,ADM_kall,ADM_gall_1d,ADM_lall   )
    real(RP) preg0_ij    (ADM_gall_1d_in,ADM_kall,ADM_gall_1d,ADM_lall   )
    real(RP) rhog0_ij    (ADM_gall_1d_in,ADM_kall,ADM_gall_1d,ADM_lall   )
    real(RP) Sr_ij       (ADM_gall_1d_in,ADM_kall,ADM_gall_1d,ADM_lall   )
    real(RP) Sw_ij       (ADM_gall_1d_in,ADM_kall,ADM_gall_1d,ADM_lall   )
    real(RP) Sp_ij       (ADM_gall_1d_in,ADM_kall,ADM_gall_1d,ADM_lall   )
#ifdef OUTPUTCHECK
    real(8) :: tmp_rhogw = 0.0_8
    real(8) :: ref_rhogw = 9527677.540549790_8
#endif
 
#ifdef USE_TIMER
    real*8 t_all_0, t_all, t_kernel_0, t_kernel
    t_all_0 = 0.0
    t_all = 0.0
    t_kernel_0 = 0.0
    t_kernel = 0.0
#endif
    PROF_INIT
    PROF_START_ALL
#ifdef USE_TIMER
    call gettod(t_all_0)
#endif

    dt = 1.0
    do l=1,ADM_lall
      do k=1,ADM_kall
        do j=1,ADM_gall_1d
          do i=1,ADM_gall_1d_in
            rhogw_ij(i,k,j,l)   = 0.0
            rhogw0_ij(i,k,j,l) = 1.0 / (i * j) / k / l / 1000
            preg0_ij(i,k,j,l) = 1.0 / (i * j) / k / l / 1000
            rhog0_ij(i,k,j,l) = 1.0  / (i * j) / k / l / 1000
            Sr_ij(i,k,j,l) = 1.0 / (i * j) / k / l / 1000
            Sw_ij(i,k,j,l) = 1.0 / (i * j) / k / l / 1000
            Sp_ij(i,k,j,l) = 1.0 / (i * j) / k / l / 1000

            VMTR_RGSGAM2_ij(i,k,j,l)= 1.0 / i / j / k / l / 1000
            VMTR_RGSGAM2H_ij(i,k,j,l)=1.0 / i / j / k / l / 1000
            VMTR_RGAMH_ij(i,k,j,l)= 1.0 / i / j / k / l / 1000
            VMTR_RGAM_ij(i,k,j,l)=1.0 / i / j / k / l / 1000
            VMTR_GSGAM2H_ij(i,k,j,l)=1.0 / i / j / k / l / 1000
            Mu_ij(i,k,j,l)=3.0
            Mc_ij(i,k,j,l)=2.0
            ML_ij(i,k,j,l)=1.0

          end do
        end do
      end do
    end do
    do k=1,ADM_kall
      GRD_rdgzh(k) = 1.0 
      GRD_afac(k) = 1.0
      GRD_bfac(k) = 1.0
    end do

    call vi_rhow_solver( &
       rhogw_ij,  rhogw_pl,  & 
       rhogw0_ij, rhogw0_pl, & 
       preg0_ij,  preg0_pl,  & 
       rhog0_ij,  rhog0_pl,  & 
       Sr_ij,     Sr_pl,     & 
       Sw_ij,     Sw_pl,     & 
       Sp_ij,     Sp_pl,     & 
       dt                 )

    PROF_START("vi_rhow_solver")
#ifdef USE_TIMER
    call gettod(t_kernel_0)
#endif
!    call timer_sta(1)
!    do i=1,110 !org
    do i=1,1
    call vi_rhow_solver( &
       rhogw_ij,  rhogw_pl,  & 
       rhogw0_ij, rhogw0_pl, & 
       preg0_ij,  preg0_pl,  & 
       rhog0_ij,  rhog0_pl,  & 
       Sr_ij,     Sr_pl,     & 
       Sw_ij,     Sw_pl,     & 
       Sp_ij,     Sp_pl,     & 
       dt                 )
    end do
!    call timer_end(1)
#ifdef USE_TIMER
    call gettod(t_kernel)
    t_kernel = t_kernel - t_kernel_0
#endif
    PROF_STOP("vi_rhow_solver")
#ifdef USE_TIMER
    call gettod(t_all)
    t_all = t_all - t_all_0    
#endif
    PROF_STOP_ALL
    PROF_FINALIZE
#ifdef USE_TIMER
    write(*,'(A, F10.7)')"WALLTIME: ",t_kernel
    write(*,'(A, F10.7)')"TOTAL TIME:",t_all
#endif

#ifdef OUTPUTCHECK
!ocl serial
    do l=1, ADM_lall
     do j =1, ADM_gall_1d
      do k=1, ADM_kall
       do i =1, ADM_gall_1d_in
       tmp_rhogw  = tmp_rhogw + dble(rhogw_ij(i,k,j,l) ** 2)
       enddo
      enddo
     enddo
    enddo
!    write(*,*) "OUTPUT rhogw=",sqrt(tmp_rhogw)
    tmp_rhogw = sqrt(tmp_rhogw)
    call report_validation(tmp_rhogw, ref_rhogw, 0.003_8)

#endif

    return
  end subroutine dynamics_step

  subroutine gettod(wctime)
    use iso_fortran_env, only: int64
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp) :: wctime
    integer(int64) :: r, c
    call system_clock(c, r)
    wctime = real(c, dp) / r
  end subroutine gettod

