#!/bin/bash
#PJM -g ra000012
#PJM -x PJM_LLIO_GFSCACHE=/vol0004
#PJM -N ML4T-fs-nicdiv-21
#PJM -L rscgrp=small
#PJM -L elapse=8:00:00
#PJM -L node=101
#PJM --mpi "max-proc-per-node=1"
# #PJM --llio localtmp-size=40Gi
#PJM -j -S

set -e

export TMPDIR=/worktmp
export LD_PRELOAD=/usr/lib/FJSVtcs/ple/lib64/libpmix.so 

module load LLVM/llvmorg-21.1.0

mpirun -n 1 python -u app.py --n-threads 2 --method FugakuEvoTADASHI --population-size 100 --max-gen 50
