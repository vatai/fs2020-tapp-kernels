/* ----- adventure region0 ---- */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef	__ARM_FEATURE_SVE
#include <arm_sve.h>
#endif

#ifndef NO_OPEN_MP
#include <omp.h>
#endif

//#include <fj_tool/fipp.h>
//#include <fjcoll.h>
#include "profiler.h"
#include "report.h"

//#include "measure.h"
//#include "matvec.h"

//#include <mpi.h>

#if 1
//  SIZE = 6*MAX_DIAG_BLOCKS
#define SIZE 5040
#endif

typedef struct {
  int nDofs;
  int nBlocks;

  int instanceId;
} MatVec;

static double vector[SIZE];
static double result[SIZE];

#define N_THREADS 12

//#define MAX_INSTANCES     32
#define MAX_INSTANCES     1
#define MAX_DIAG_BLOCKS  840
#define MAX_OFF_DIAG_BLOCKS ((MAX_DIAG_BLOCKS * (MAX_DIAG_BLOCKS - 1)) / 2)

static double DiagComponents[MAX_INSTANCES][MAX_DIAG_BLOCKS][6][6];
static double OffDiagComponents[MAX_INSTANCES][MAX_OFF_DIAG_BLOCKS][6][6];

static MatVec TheMatVec;

void MatVec_clearMatrix
(MatVec *self)
{
  int instanceId = self->instanceId;
  int nBlocks = self->nBlocks;
  int nOffDiagBlocks;
  int iBlock;
  /* add */
  int jBlock;
  /* add end. */
  int i, j;

#pragma omp parallel
#pragma omp for private(iBlock,i,j)
  for (iBlock = 0;
       iBlock <= nBlocks;
       iBlock++) {

/* original.
#pragma omp parallel
#pragma omp for private(j)
end. */
#pragma loop xfill
    for (i = 0; i < 6; i++) {
      for (j = 0; j < 6; j++) {


        /* original.
        DiagComponents[instanceId][iBlock][i][j] = 0.0;
        end. */
        //DiagComponents[instanceId][iBlock][i][j] = 1.0e0;
        /* alternative for debug. */
        DiagComponents[instanceId][iBlock][i][j] = (double)(iBlock*6 + j + 1);
      }
    }
  }

  nOffDiagBlocks = (nBlocks * (nBlocks - 1)) / 2;
  /* original.
  for (iBlock = 0;
       iBlock < nOffDiagBlocks;
       iBlock++) {

#pragma loop xfill
    for (i = 0; i < 6; i++) {
      for (j = 0; j < 6; j++) {
        OffDiagComponents[instanceId][iBlock][i][j] = 0.0;
      }
    }
  }
  end. */
#pragma omp parallel
#pragma omp for private(iBlock,jBlock,i,j)
  for (iBlock = 0; iBlock < nBlocks; iBlock++) {
    int offset = (iBlock * (iBlock - 1)) / 2;
    for (jBlock = 0; jBlock < iBlock; jBlock++) {
      for (i = 0; i < 6; i++) {
        for (j = 0; j < 6; j++) {
          OffDiagComponents[instanceId][offset + jBlock][i][j] =
          //  1.0e0;
          /* alternative for debug. */
            (double)(jBlock*6 + j + 1);
        }
      }
    } /* jBlock */
  } /* iBlock */

}
void MatVec_product
(MatVec *self,
 double vector[],
 double result_OUT[])
{
  //padding.  static double tmp[N_THREADS][MAX_DIAG_BLOCKS][6];
  static double tmp[N_THREADS][MAX_DIAG_BLOCKS][8];
  static double vector2[MAX_DIAG_BLOCKS][6];
  //padding.  static double result2[MAX_DIAG_BLOCKS][6];
  static double result2[MAX_DIAG_BLOCKS][8];
  /* for SIMD-loop */
  //padding.  static double vec_m[6][6];
  static double vec_m[6][9];
  /* end */
  int instanceId = self->instanceId;
  int nDofs = self->nDofs;
  //int nBlocks = self->nBlocks;
  int nBlocks = MAX_DIAG_BLOCKS;
  int nFullBlocks = nDofs / 6;
  int iBlock;
  int iDof;
  int iThread;

#pragma omp parallel
#pragma omp  for
  for (iThread = 0; iThread < N_THREADS; iThread++) {
    int iBlock;

    for (iBlock = 0; iBlock < nBlocks; iBlock++) {
      tmp[iThread][iBlock][0] = 0.0;
      tmp[iThread][iBlock][1] = 0.0;
      tmp[iThread][iBlock][2] = 0.0;
      tmp[iThread][iBlock][3] = 0.0;
      tmp[iThread][iBlock][4] = 0.0;
      tmp[iThread][iBlock][5] = 0.0;
    }
  }

#pragma omp parallel
#pragma omp  for
  for (iBlock = 0; iBlock < nFullBlocks; iBlock++) {
    vector2[iBlock][0] = vector[iBlock * 6 + 0];
    vector2[iBlock][1] = vector[iBlock * 6 + 1];
    vector2[iBlock][2] = vector[iBlock * 6 + 2];
    vector2[iBlock][3] = vector[iBlock * 6 + 3];
    vector2[iBlock][4] = vector[iBlock * 6 + 4];
    vector2[iBlock][5] = vector[iBlock * 6 + 5];
  }

  if (nFullBlocks * 6 < nDofs) {

    for (iDof = 0;
         iDof < nDofs % 6;
         iDof++) {
      vector2[nFullBlocks][iDof] = vector[nFullBlocks * 6 + iDof];
    }

    for (iDof = nDofs % 6;
         iDof < 6;
         iDof++) {
      vector2[nFullBlocks][iDof] = 0.0;
    }
  }

  int threadId = 0;
  int Chunk = 12;

  int ncycle = (int)( (int)( (nBlocks - 1)/Chunk )/N_THREADS ) + 1;
#pragma omp parallel
#pragma omp for private(vec_m)
  for (threadId = 0 ; threadId < N_THREADS ; threadId++ ){

  int i, ic;
  int iasc = threadId + 1;
  int idsc = N_THREADS - threadId;

  for(ic = 0; ic < ncycle; ic++){
    int ithr = iasc;
    if(ic%2 == 1){ ithr = idsc; }
    int ioff = nBlocks - (N_THREADS*ic + ithr)*Chunk;    /* descending order */

  for (i = 0 ; i < Chunk; i++) {
    int iBlock = ioff + i;
    if(iBlock < 0) continue;                             /* descending order */

    int offset = (iBlock * (iBlock - 1)) / 2;
    int iv;
    int jv;

    for ( jv = 0; jv < 6; jv++ ){ for ( iv = 0; iv < 6; iv++ ){ vec_m[jv][iv] = 0.0e0; } }

    double v_i0 = vector2[iBlock][0];
    double v_i1 = vector2[iBlock][1];
    double v_i2 = vector2[iBlock][2];
    double v_i3 = vector2[iBlock][3];
    double v_i4 = vector2[iBlock][4];
    double v_i5 = vector2[iBlock][5];

    int jBlock;

#ifdef	__ARM_FEATURE_SVE
    svbool_t p0 = svptrue_pat_b64(SV_VL6);

    svfloat64_t v_i0_v = svdup_n_f64(v_i0);
    svfloat64_t v_i1_v = svdup_n_f64(v_i1);
    svfloat64_t v_i2_v = svdup_n_f64(v_i2);
    svfloat64_t v_i3_v = svdup_n_f64(v_i3);
    svfloat64_t v_i4_v = svdup_n_f64(v_i4);
    svfloat64_t v_i5_v = svdup_n_f64(v_i5);

    svfloat64_t vec_m_0 = svld1_f64(p0, &vec_m[0][0]);
    svfloat64_t vec_m_1 = svld1_f64(p0, &vec_m[1][0]);
    svfloat64_t vec_m_2 = svld1_f64(p0, &vec_m[2][0]);
    svfloat64_t vec_m_3 = svld1_f64(p0, &vec_m[3][0]);
    svfloat64_t vec_m_4 = svld1_f64(p0, &vec_m[4][0]);
    svfloat64_t vec_m_5 = svld1_f64(p0, &vec_m[5][0]);

    for (jBlock = 0; jBlock < iBlock; jBlock++) {
	/* these prefetch intrinsics are local to the development vehicle, and not available for general use.
	 * Alternatively, use prefetch directives available on the tested systems
      __builtin_fj_prefetch(&OffDiagComponents[instanceId][offset + jBlock][0][0]+1024, 0, 1, 1);
      __builtin_fj_prefetch(&OffDiagComponents[instanceId][offset + jBlock][0][0]+2048, 0, 2, 1);
	*/

      /* lower-triangular part: row-wised SIMD */
      svfloat64_t v_j    = svld1_f64(p0, &vector2[jBlock][0]);
      svfloat64_t m_i0_j = svld1_f64(p0, &OffDiagComponents[instanceId][offset + jBlock][0][0]);
      svfloat64_t m_i1_j = svld1_f64(p0, &OffDiagComponents[instanceId][offset + jBlock][1][0]);
      svfloat64_t m_i2_j = svld1_f64(p0, &OffDiagComponents[instanceId][offset + jBlock][2][0]);
      svfloat64_t m_i3_j = svld1_f64(p0, &OffDiagComponents[instanceId][offset + jBlock][3][0]);
      svfloat64_t m_i4_j = svld1_f64(p0, &OffDiagComponents[instanceId][offset + jBlock][4][0]);
      svfloat64_t m_i5_j = svld1_f64(p0, &OffDiagComponents[instanceId][offset + jBlock][5][0]);

        //vec_m[0][iv] += m_i0_j * v_j;
        //vec_m[1][iv] += m_i1_j * v_j;
        //vec_m[2][iv] += m_i2_j * v_j;
        //vec_m[3][iv] += m_i3_j * v_j;
        //vec_m[4][iv] += m_i4_j * v_j;
        //vec_m[5][iv] += m_i5_j * v_j;
      vec_m_0 = svmla_f64_m(p0, vec_m_0, m_i0_j, v_j);
      vec_m_1 = svmla_f64_m(p0, vec_m_1, m_i1_j, v_j);
      vec_m_2 = svmla_f64_m(p0, vec_m_2, m_i2_j, v_j);
      vec_m_3 = svmla_f64_m(p0, vec_m_3, m_i3_j, v_j);
      vec_m_4 = svmla_f64_m(p0, vec_m_4, m_i4_j, v_j);
      vec_m_5 = svmla_f64_m(p0, vec_m_5, m_i5_j, v_j);

      /* upper-triangular part: transposed column-wised SIMD */
        //tmp[threadId][jBlock][iv]
        //          += m_i0_j * v_i0 + m_i1_j * v_i1 + m_i2_j * v_i2
        //           + m_i3_j * v_i3 + m_i4_j * v_i4 + m_i5_j * v_i5;

      svfloat64_t tmp_v = svld1_f64(p0, &tmp[threadId][jBlock][0]);
      tmp_v = svmla_f64_m(p0, tmp_v, m_i0_j, v_i0_v);
      tmp_v = svmla_f64_m(p0, tmp_v, m_i1_j, v_i1_v);
      tmp_v = svmla_f64_m(p0, tmp_v, m_i2_j, v_i2_v);
      tmp_v = svmla_f64_m(p0, tmp_v, m_i3_j, v_i3_v);
      tmp_v = svmla_f64_m(p0, tmp_v, m_i4_j, v_i4_v);
      tmp_v = svmla_f64_m(p0, tmp_v, m_i5_j, v_i5_v);
      svst1_f64(p0, &tmp[threadId][jBlock][0], tmp_v);
    }    /* jblock */

    svst1_f64(p0, &vec_m[0][0], vec_m_0);
    svst1_f64(p0, &vec_m[1][0], vec_m_1);
    svst1_f64(p0, &vec_m[2][0], vec_m_2);
    svst1_f64(p0, &vec_m[3][0], vec_m_3);
    svst1_f64(p0, &vec_m[4][0], vec_m_4);
    svst1_f64(p0, &vec_m[5][0], vec_m_5);
#else

    for (jBlock = 0; jBlock < iBlock; jBlock++) {
	/* these prefetch intrinsics are local to the development vehicle, and not available for general use.
	 * Alternatively, use prefetch directives available on the tested systems
      __builtin_fj_prefetch(&OffDiagComponents[instanceId][offset + jBlock][0][0]+1024, 0, 1, 1);
      __builtin_fj_prefetch(&OffDiagComponents[instanceId][offset + jBlock][0][0]+2048, 0, 2, 1);
	*/

      /* lower-triangular part: row-wised SIMD */
#pragma loop simd_redundant_vl 6
      for ( iv = 0; iv < 6; iv++ ){
        double v_j = vector2[jBlock][iv];
        double m_i0_j = OffDiagComponents[instanceId][offset + jBlock][0][iv];
        double m_i1_j = OffDiagComponents[instanceId][offset + jBlock][1][iv];
        double m_i2_j = OffDiagComponents[instanceId][offset + jBlock][2][iv];
        double m_i3_j = OffDiagComponents[instanceId][offset + jBlock][3][iv];
        double m_i4_j = OffDiagComponents[instanceId][offset + jBlock][4][iv];
        double m_i5_j = OffDiagComponents[instanceId][offset + jBlock][5][iv];
        vec_m[0][iv] += m_i0_j * v_j;
        vec_m[1][iv] += m_i1_j * v_j;
        vec_m[2][iv] += m_i2_j * v_j;
        vec_m[3][iv] += m_i3_j * v_j;
        vec_m[4][iv] += m_i4_j * v_j;
        vec_m[5][iv] += m_i5_j * v_j;
      }
      /* upper-triangular part: transposed column-wised SIMD */
#pragma loop simd_redundant_vl 6
      for ( iv = 0; iv < 6; iv++ ){
        double m_i0_j = OffDiagComponents[instanceId][offset + jBlock][0][iv];
        double m_i1_j = OffDiagComponents[instanceId][offset + jBlock][1][iv];
        double m_i2_j = OffDiagComponents[instanceId][offset + jBlock][2][iv];
        double m_i3_j = OffDiagComponents[instanceId][offset + jBlock][3][iv];
        double m_i4_j = OffDiagComponents[instanceId][offset + jBlock][4][iv];
        double m_i5_j = OffDiagComponents[instanceId][offset + jBlock][5][iv];
        tmp[threadId][jBlock][iv]
                  += m_i0_j * v_i0 + m_i1_j * v_i1 + m_i2_j * v_i2
                   + m_i3_j * v_i3 + m_i4_j * v_i4 + m_i5_j * v_i5;
      }

    }    /* jblock */
#endif

    {    /* diagonal block */
      double m_i0_j0 = DiagComponents[instanceId][iBlock][0][0];

      double m_i1_j0 = DiagComponents[instanceId][iBlock][1][0];
      double m_i1_j1 = DiagComponents[instanceId][iBlock][1][1];

      double m_i2_j0 = DiagComponents[instanceId][iBlock][2][0];
      double m_i2_j1 = DiagComponents[instanceId][iBlock][2][1];
      double m_i2_j2 = DiagComponents[instanceId][iBlock][2][2];

      double m_i3_j0 = DiagComponents[instanceId][iBlock][3][0];
      double m_i3_j1 = DiagComponents[instanceId][iBlock][3][1];
      double m_i3_j2 = DiagComponents[instanceId][iBlock][3][2];
      double m_i3_j3 = DiagComponents[instanceId][iBlock][3][3];

      double m_i4_j0 = DiagComponents[instanceId][iBlock][4][0];
      double m_i4_j1 = DiagComponents[instanceId][iBlock][4][1];
      double m_i4_j2 = DiagComponents[instanceId][iBlock][4][2];
      double m_i4_j3 = DiagComponents[instanceId][iBlock][4][3];
      double m_i4_j4 = DiagComponents[instanceId][iBlock][4][4];

      double m_i5_j0 = DiagComponents[instanceId][iBlock][5][0];
      double m_i5_j1 = DiagComponents[instanceId][iBlock][5][1];
      double m_i5_j2 = DiagComponents[instanceId][iBlock][5][2];
      double m_i5_j3 = DiagComponents[instanceId][iBlock][5][3];
      double m_i5_j4 = DiagComponents[instanceId][iBlock][5][4];
      double m_i5_j5 = DiagComponents[instanceId][iBlock][5][5];

      double m_i0_j1 = m_i1_j0;

      double m_i0_j2 = m_i2_j0;
      double m_i1_j2 = m_i2_j1;

      double m_i0_j3 = m_i3_j0;
      double m_i1_j3 = m_i3_j1;
      double m_i2_j3 = m_i3_j2;

      double m_i0_j4 = m_i4_j0;
      double m_i1_j4 = m_i4_j1;
      double m_i2_j4 = m_i4_j2;
      double m_i3_j4 = m_i4_j3;

      double m_i0_j5 = m_i5_j0;
      double m_i1_j5 = m_i5_j1;
      double m_i2_j5 = m_i5_j2;
      double m_i3_j5 = m_i5_j3;
      double m_i4_j5 = m_i5_j4;

      double v_i0 = vector2[iBlock][0];
      double v_i1 = vector2[iBlock][1];
      double v_i2 = vector2[iBlock][2];
      double v_i3 = vector2[iBlock][3];
      double v_i4 = vector2[iBlock][4];
      double v_i5 = vector2[iBlock][5];
      vec_m[0][0] += m_i0_j0 * v_i0 + m_i0_j1 * v_i1 + m_i0_j2 * v_i2
                  + m_i0_j3 * v_i3 + m_i0_j4 * v_i4 + m_i0_j5 * v_i5;
      vec_m[1][0] += m_i1_j0 * v_i0 + m_i1_j1 * v_i1 + m_i1_j2 * v_i2
                  + m_i1_j3 * v_i3 + m_i1_j4 * v_i4 + m_i1_j5 * v_i5;
      vec_m[2][0] += m_i2_j0 * v_i0 + m_i2_j1 * v_i1 + m_i2_j2 * v_i2
                  + m_i2_j3 * v_i3 + m_i2_j4 * v_i4 + m_i2_j5 * v_i5;
      vec_m[3][0] += m_i3_j0 * v_i0 + m_i3_j1 * v_i1 + m_i3_j2 * v_i2
                  + m_i3_j3 * v_i3 + m_i3_j4 * v_i4 + m_i3_j5 * v_i5;
      vec_m[4][0] += m_i4_j0 * v_i0 + m_i4_j1 * v_i1 + m_i4_j2 * v_i2
                  + m_i4_j3 * v_i3 + m_i4_j4 * v_i4 + m_i4_j5 * v_i5;
      vec_m[5][0] += m_i5_j0 * v_i0 + m_i5_j1 * v_i1 + m_i5_j2 * v_i2
                  + m_i5_j3 * v_i3 + m_i5_j4 * v_i4 + m_i5_j5 * v_i5;

    }    /* diagonal block */

#pragma loop simd_redundant_vl 5
    for ( iv = 1; iv < 6; iv++ ){
      vec_m[0][0] += vec_m[0][iv];
      vec_m[1][0] += vec_m[1][iv];
      vec_m[2][0] += vec_m[2][iv];
      vec_m[3][0] += vec_m[3][iv];
      vec_m[4][0] += vec_m[4][iv];
      vec_m[5][0] += vec_m[5][iv];
    }

    for ( iv = 0; iv < 6; iv++ ){
      result2[iBlock][iv] = vec_m[iv][0];
    }

  }  /* i : within chunk of row blocks */
  }  /* ic : reciprocal cyclic assignment of a chunk of blocks for threads */
  }  /* threadId : threads */


  int div = (nBlocks - 1) / N_THREADS + 1;
  int ithr;
#pragma omp parallel 
#pragma omp for private(iThread)
  for (ithr = 0; ithr < N_THREADS; ithr++) {
    int iBend = (ithr + 1) * div;
    if (iBend > nBlocks){ iBend = nBlocks; }
    for (iThread = 0; iThread < N_THREADS; iThread++) {
    int iBlock;
    int iv;
    for (iBlock = ithr*div; iBlock < iBend; iBlock++) {
      for (iv = 0; iv <6; iv++) {
        result2[iBlock][iv] += tmp[iThread][iBlock][iv];
      }
    }
    }
  }

#pragma omp parallel 
#pragma omp for
  for (iBlock = 0; iBlock < nFullBlocks; iBlock++) {
    result_OUT[iBlock * 6 + 0] = result2[iBlock][0];
    result_OUT[iBlock * 6 + 1] = result2[iBlock][1];
    result_OUT[iBlock * 6 + 2] = result2[iBlock][2];
    result_OUT[iBlock * 6 + 3] = result2[iBlock][3];
    result_OUT[iBlock * 6 + 4] = result2[iBlock][4];
    result_OUT[iBlock * 6 + 5] = result2[iBlock][5];
  }

  if (nFullBlocks * 6 < nDofs) {

    for (iDof = 0;
         iDof < nDofs % 6;
         iDof++) {
      result_OUT[nFullBlocks * 6 + iDof] = result2[nFullBlocks][iDof];
    }
  }

#ifdef MEASURE
  Measure_AddOp (M_ADD_OP, nDofs * nDofs);
  Measure_AddOp (M_MUL_OP, nDofs * nDofs);
#endif
}

int main()
{

  int iIteration;
  int nDofs = SIZE;
  int i;
  /* add */
  int j;
  /* add end */
  PROF_INIT;
  PROF_START_ALL;
  TheMatVec.instanceId = 0;
  TheMatVec.nDofs = SIZE;
  TheMatVec.nBlocks = MAX_DIAG_BLOCKS;

  MatVec_clearMatrix(&TheMatVec);
  for(i=0;i<nDofs;i++) {
     vector[i] = 1.0e0;
     //vector[i] = (double)i + 1.0e0;
  }
  for (iIteration = 0; iIteration < 2; iIteration++) {
  if(iIteration !=0 )  {
      PROF_START("region0");
  }
  /*
  start_collection("region0_tune");
  */
  //for (iIteration = 0; iIteration < 100; iIteration++) {
  //for (iIteration = 0; iIteration < 4; iIteration++) {
  //for (iIteration = 0; iIteration < 1; iIteration++) {
  /* iteration for simulation. */
    MatVec_product
      (&TheMatVec,
       vector,
       result);
  //}
  /*
  stop_collection("region0_tune");
  */
  if(iIteration !=0 )  {
      PROF_STOP("region0");
  }
}
  PROF_STOP_ALL;
  PROF_FINALIZE;
  /* add. */
  double sum = 0.0e0;
  for(i=0; i<SIZE; i++){
    sum += result[i];
  }
  //printf("\n");
  //printf(" Result Validation. Value: %0.10f",sum);
  //printf("\n");

  report_validation(sum,4.268738964000000e+10,1.0e-13);

  /* add for debug. */
  /*
  for(i=0; i<SIZE; i=i+8){
    int jmax = i+8;
    if(jmax > SIZE) { jmax = SIZE; }
    for(j=i; j<jmax; j++){
      printf(" %0.10f",result[j]);
    }
    printf("\n");
  }
  */
  /* add end */
}
/* ----- end of adventure region0 ---- */
