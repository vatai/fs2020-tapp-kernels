#!/bin/bash
#PJM -N QCD-JINV
#PJM --rsc-list "rscgrp=small"
#PJM --rsc-list "elapse=00:10:00"
#PJM --mpi "max-proc-per-node=1"
#PJM --rsc-list "node=1"
#PJM -j
#PJM -S

# module list
# date
# hostname
# 
# TMPDIR=${HOME}/tmp/check_QCD-JINV
# mkdir -p ${TMPDIR}
# cd ${TMPDIR}/
# if [ $? != 0 ] ; then echo '@@@ Directory error @@@'; exit; fi
# pwd
# rm *
# 
# SRC_DIR=${HOME}/fs2020_kernels/src/QCD.jinv_ddd_in_s_.20170720
# cp -rp ${SRC_DIR}/* ./

OPTIMIZE="-fopenmp"
DEFS="-DRDC -DVLENS=16 -DEML_LIB -DPREFETCH -D_CHECK_SIM -DTARGET_JINV"
CXXFLAGS="${OPTIMIZE} ${DEFS} -std=gnu++11 -Icommon/include "
CFLAGS="${OPTIMIZE} ${DEFS} -std=c99 -Icommon/include -DDISABLE_VALIDATION "

set -x

for i in \
bicgstab_precdd_s.cc clover_s.cc ddd_in_s.cc ddd_out_s.cc main.cc qws.cc static_solver.cc
do
clang++ -c ${CXXFLAGS} ${i}
done
clang -c ${CFLAGS} tools.c
clang -c ${CFLAGS} common/src/report.c 

clang++ ${CXXFLAGS} *.o

export OMP_NUM_THREADS=12
time ./a.out
exit
