#include "profiler.h"
program main
  !!$ use omp_lib
  implicit none
  external omp_get_thread_num
  integer omp_get_thread_num
  integer(kind=  4),parameter :: N = 393216
!  integer(kind=  4),parameter :: M = 12
  integer(kind=  4),parameter :: M = 24
  integer(kind=  4),parameter :: RKP = 8
  integer(kind=4),parameter :: BLEN = 10000
!  real   (kind=RKP),dimension(:,:),allocatable :: matrix_H
!  real   (kind=RKP),dimension(:,:),allocatable :: matrix_L
!  real   (kind=RKP),dimension(:  ),allocatable :: vector_i
!  real   (kind=RKP),dimension(:  ),allocatable :: vector_o
  real   (kind=RKP),dimension(N,M) :: matrix_H
  real   (kind=RKP),dimension(N,M) :: matrix_L
  real   (kind=RKP),dimension(N  ) :: vector_i
  real   (kind=RKP),dimension(  M) :: vector_o
  integer(kind=  4) :: itr,i,j
  integer(kind=  8) :: iaddr1, iaddr2, iaddr3

  integer(kind=4),parameter :: MAX_NUM_THREADS = 16
  integer(kind=4) :: thread_num
  real   (kind=8),dimension(BLEN+1, 0:MAX_NUM_THREADS-1):: t1_H,t1_L
  real   (kind=8) :: t1_H_ss,t1_L_ss
  !
!  allocate(matrix_H(N,M))
!  allocate(matrix_L(N,M))
!  allocate(vector_i(N  ))
!  allocate(vector_o(  M))
  !
  PROF_INIT
  PROF_START_ALL

  do j=1,M
     do i=1,N
        if(i.eq.1) then
           matrix_H(i,j)=real( 1.0E+00_8,kind=RKP)
           matrix_L(i,j)=real( 0.0E+00_8,kind=RKP)
        else if(i.eq.N) then
           matrix_H(i,j)=real(-1.0E+00_8,kind=RKP)
           matrix_L(i,j)=real( 0.0E+00_8,kind=RKP)
        else
           matrix_H(i,j)=real( 1.0E-9_8,kind=RKP)
           matrix_L(i,j)=real( 1.01E-17_8,kind=RKP)
           matrix_H(i,j) = matrix_H(i,j) * (mod(i,5)+1)
           matrix_L(i,j) = matrix_L(i,j) * (mod(i,5)+1)
        end if
     end do
  end do
  do i=1,N
     vector_i(i)=real(1.0E+00_8,kind=RKP)
     vector_i(i) = vector_i(i)  * (mod(i,5)+1)
  end do
  t1_H = 0.0_8
  t1_L = 0.0_8
  !
!  do itr=1,1000
  do itr=1,2
      if(itr > 1) then
          PROF_START("region1")
      endif
     !$omp parallel default(shared), private(j, thread_num)
      thread_num = omp_get_thread_num()
     !$omp do
     do j=1,M
!        call DD_InnerProduct(N,matrix_H(:,j),matrix_L(:,j),vector_i(:),vector_o(j))
        call DD_InnerProduct(N,matrix_H(:,j),matrix_L(:,j),vector_i(:),vector_o(j),BLEN,t1_H(:,thread_num),t1_L(:,thread_num))
     end do
     !$omp end do
     !$omp end parallel
      if(itr > 1) then
          PROF_STOP("region1")
      endif
  end do

  call get_ss_r8(t1_H, (BLEN+1)*MAX_NUM_THREADS, t1_H_ss)
  call get_ss_r8(t1_L, (BLEN+1)*MAX_NUM_THREADS, t1_L_ss)

  ! the reference values are valid only when #threads = 12
  call report_validation(t1_H_ss, 4.800000000000240e+01_8, 0_8)
  call report_validation(t1_L_ss, 2.478474534255483e-44_8, 0_8)
  
  !
!  write(6,'(e30.20)') real(vector_o(1),kind=8)
  !
  iaddr1=loc(matrix_H(1,1))
  iaddr2=loc(matrix_L(1,1))
  iaddr3=loc(vector_i(1)  )
  !
!  call iaddr_print(iaddr1)
!  call iaddr_print(iaddr2)
!  call iaddr_print(iaddr3)
  !
!  deallocate(matrix_H,matrix_L)
!  deallocate(vector_i,vector_o)
  !
  PROF_STOP_ALL
  PROF_FINALIZE
end program main
!
!
!subroutine DD_InnerProduct(N,matrix_H,matrix_L,vector_i,vector_o)
subroutine DD_InnerProduct(N,matrix_H,matrix_L,vector_i,vector_o,BLEN,t1_H,t1_L)
  implicit none
  integer(kind=4),             intent(IN)  :: N
  real   (kind=8),dimension(N),intent(IN)  :: matrix_H, matrix_L
  real   (kind=8),dimension(N),intent(IN)  :: vector_i
  real   (kind=8)             ,intent(OUT) :: vector_o
  integer(kind=4),intent(IN) :: BLEN
  !
!  integer(kind=4),parameter :: BLEN = 10000
  real   (kind=8),parameter :: QD_SPLITTER = 134217729.0E0_8
 ! real   (kind=8),dimension(:),allocatable :: t1_H,t1_L,t2_H,t2_L
  real   (kind=8),dimension(BLEN    +1):: t1_H,t1_L
  real   (kind=8),dimension(BLEN/2+1+1):: t2_H,t2_L
  real   (kind=8) :: vector_H,vector_L
  real   (kind=8) :: p3,temp0
  real   (kind=8) :: aa_hi,aa_lo
  real   (kind=8) :: temp1
  real   (kind=8) :: bb_hi,bb_lo
  real   (kind=8) :: p1_L1,p1_L2
  real   (kind=8) :: p2_L1,p2_L2,p2_L3
  real   (kind=8) :: a0,a1
  real   (kind=8) :: b0,b1,bb0
  real   (kind=8) :: e_L1,e_L2,e_L3
  real   (kind=8) :: s0,s1,s_L1,s_L2
  integer(kind=4) :: h,i,k,ii,nn
  !
!  allocate(t1_H(BLEN    +1),t1_L(BLEN    +1))
!  allocate(t2_H(BLEN/2+1+1),t2_L(BLEN/2+1+1))
  vector_H = 0.0E0_8
  vector_L = 0.0E0_8
!ocl unroll(3)
  do ii=1,N,BLEN
     nn=min(BLEN,N-ii+1)
     do i=1,nn
        a0 = matrix_H(ii+i-1)
        a1 = matrix_L(ii+i-1)
        b0 = vector_i(ii+i-1)
        ! b1 = 0.0E0_8
        p3 = a0 * b0
        temp0 = QD_SPLITTER * a0
        aa_hi = temp0 - (temp0 - a0)
        aa_lo = a0 - aa_hi
        temp1 = QD_SPLITTER * b0
        bb_hi = temp1 - (temp1 - b0)
        bb_lo = b0 - bb_hi
        p2_L1 = ((aa_hi * bb_hi - p3) + aa_hi * bb_lo + aa_lo * bb_hi) + aa_lo * bb_lo
        p1_L1 = p3
        p2_L2 = p2_L1 + (a1 * b0)  ! <-- p2_L2 = p2_L1 + (a0 * b1 + a1 * b0)
        s1 = p1_L1 + p2_L2
        p2_L3 = p2_L2 - (s1 - p1_L1)
        p1_L2 = s1
        t1_H(i) = p1_L2
        t1_L(i) = p2_L3
     end do
     !
  end do
  vector_o = vector_H
!  deallocate(t1_H,t1_L,t2_H,t2_L)
  !
  return
end subroutine DD_InnerProduct


subroutine get_ss_r8(array, length, sum_square)
    real(kind=8),dimension(length) :: array
    real(kind=8) :: sum_square
    integer :: length

    sum_square = 0.0
    do i=1,length
    sum_square = sum_square + array(i)*array(i)
    end do

end subroutine get_ss_r8

subroutine report_validation (computed_result, expected_result, error_bar)
    real(kind=8) :: computed_result, expected_result, error_bar
    real(kind=8) :: error_norm
         
    error_norm = computed_result/expected_result

    if ( 0.999 < error_norm .and. error_norm < 1.001 ) then
        write(*,*) "the computed result seems to be OK."
    else
     write(*,*) "the computed result is not close enough to the expected value."     
    endif
      
end subroutine


