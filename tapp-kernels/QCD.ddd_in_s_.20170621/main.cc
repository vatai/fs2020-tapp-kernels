/*
!===============================================================================
!
! Copyright (C) 2014,2015 Yoshifumi Nakamura
!
! _CODENAME_ = ???
!
! This file is part of _CODENAME_
!
! _CODENAME_ is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! _CODENAME_ is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with _CODENAME_.  If not, see <http://www.gnu.org/licenses/>.
!
!-----------------------------------------------------------------------------
*/
#include "qws.h"
#include "qwsintrin.h"
#include "report.h"
#include "profiler.h"

#include <stdio.h>
#include <ctime>
#include <time.h>

extern int vold, vols, rank, nx, ny, nz, nt, nxh, nxd, nxs;
extern double kappa2, kappa, mkappa;
extern "C" void qws_init_(int* lx,  int* ly, int* lz, int* lt, int* npe_f, int* fbc_f, int* pce_f, int* pco_f, int* block_size);
extern "C" void bicgstab_dd_mix_(scd_t* x, scd_t* b, double *tol, int* conviter, int* maxiter, double* tol_s, int* maxiter_s, int* nsap, int* nm);

float* result_arr;
int result_arr_len;

int main()
{
  PROF_INIT;
  PROF_START_ALL;
  
  int lx=NX;
  int ly=NY;
  int lz=NZ;
  int lt=NT;
  int npe_f[4];
  npe_f[0]=1;
  npe_f[1]=1;
  npe_f[2]=1;
  npe_f[3]=1;

  double tol  = -1.0;
  double tol_s= -1.0;

  int fbc_f[4];
  fbc_f[0]=1;
  fbc_f[1]=1;
  fbc_f[2]=1;
  fbc_f[3]=-1;
  int pce_f=0;
  int pco_f=1;
  
  int block_size[4];
  block_size[0]=1;
  block_size[1]=2;
  block_size[2]=2;
  block_size[3]=2;
  
  qws_init_(&lx, &ly, &lz, &lt, npe_f, fbc_f, &pce_f, &pco_f, block_size);

  kappa = (double)0.05;
  kappa2 = kappa * kappa;
  mkappa = - kappa;

  int iter;
#ifndef RDC
  int maxiter  =6;
  int maxiter_s=50;
#else
  int maxiter  =3;
  int maxiter_s=1;
#endif  
  int nsap, nm;
  nsap = 4;
  nm = 2;
  struct timespec start, end;
  clock_gettime(CLOCK_REALTIME, &start);
  bicgstab_dd_mix_(0, 0, &tol, &iter, &maxiter, &tol_s, &maxiter_s, &nsap, &nm);
  clock_gettime(CLOCK_REALTIME, &end);
  double second = (end.tv_sec - start.tv_sec) +  (end.tv_nsec - start.tv_nsec) * 1.0e-9;
  printf("WALLTIME: %fs\n", second);


#ifdef DS_TO_DOUBLE
  report_validation(get_ss_r8(result_arr, result_arr_len), -4.180054191179794e+307, 0.01);
#else
  report_validation(get_ss_r4(result_arr, result_arr_len),  4.747487304687500e+03, 0.01);
#endif

  PROF_STOP_ALL;
  PROF_FINALIZE;

  return 0;
}
