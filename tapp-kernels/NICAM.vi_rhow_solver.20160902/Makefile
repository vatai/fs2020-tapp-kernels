FC = flang
CFLAGS = -Icommon/include  -fopenmp
FCFLAGS = -Icommon/include -cpp -DUSE_TIMER -fopenmp
APP = postK_nicam_vi_rhow_solver

postK_nicam_vi_rhow_solver.o: postK_nicam_vi_rhow_solver.f90
	$(FC) $(FCFLAGS) $^ -c -o $@

main.o: main.f90
	$(FC) $(FCFLAGS) $^ -c -o $@

$(APP).x: main.o $(APP).o common/src/report.o
	$(FC) $(FCFLAGS) $^ -o $@

scop-detect:
	scop-detector . -e f -a flang $(FCFLAGS)

