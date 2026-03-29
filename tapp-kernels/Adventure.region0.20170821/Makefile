CC = clang
CFLAGS = -Icommon/include -fopenmp
APP = adventure_kernel_region0_tune4_arm_pad-acle
$(APP).x: $(APP).o common/src/report.o
	$(CC) $(CFLAGS) $^ -o $@
