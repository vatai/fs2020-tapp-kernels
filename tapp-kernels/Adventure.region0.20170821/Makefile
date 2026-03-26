CC = clang
CFLAGS = -Icommon/include  -fopenmp # --std=c99 -D__ARM_FEATURE_SVE=1 
APP = adventure_kernel_region0_tune4_arm_pad-acle
$(APP).x: $(APP).o common/src/report.o
	$(CC) $(CFLAGS) $^ -o $@

scop-detect:
	scop-detector -c . -e c $(CFLAGS)
