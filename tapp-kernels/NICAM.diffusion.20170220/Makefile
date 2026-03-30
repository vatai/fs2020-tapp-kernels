FC = flang
CFLAGS = -Icommon/include  -fopenmp
FCFLAGS = -Icommon/include  -fopenmp
APP = postK_nicam_diffusion
$(APP).x: $(APP).f90 common/src/report.o
	$(FC) $(FCFLAGS) $^ -o $@

scop-detect:
	scop-detector . -e f -a flang -Icommon/include
