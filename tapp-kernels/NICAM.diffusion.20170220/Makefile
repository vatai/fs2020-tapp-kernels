FC = flang
CFLAGS = -Icommon/include  -fopenmp
FCFLAGS = -Icommon/include  -fopenmp
APP = postK_nicam_diffusion

$(APP).x: $(APP).o common/src/report.o
	$(FC) $(FCFLAGS) $^ -o $@

postK_nicam_diffusion.o: postK_nicam_diffusion.f90
	$(FC) $(FCFLAGS) $^ -c -o $@


scop-detect:
	scop-detector . -e f -a flang -Icommon/include

