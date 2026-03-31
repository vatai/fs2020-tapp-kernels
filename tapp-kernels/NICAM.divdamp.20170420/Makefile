FC = flang
CFLAGS = -Icommon/include  -fopenmp
FCFLAGS = -Icommon/include  -fopenmp
APP = postK_nicam_divdamp_ijsplit_tune03

$(APP).x: $(APP).o common/src/report.o
	$(FC) $(FCFLAGS) $^ -o $@

postK_nicam_divdamp_ijsplit_tune03.o: postK_nicam_divdamp_ijsplit_tune03.f90
	$(FC) $(FCFLAGS) $^ -c -o $@


scop-detect:
	scop-detector . -e f -a flang $(FCFLAGS)


