CC = clang
CFLAGS = -Icommon/include
FC = flang
FCFLAGS = -Icommon/include  -fopenmp -DREAL_4=real -DMAXCOLOR=100
APP = cal_amat_tet10_s_inline_simd_color_loopdiv_r3_proposed1

cal_amat_tet10_s_inline_simd_color_loopdiv_r3_proposed1.o: cal_amat_tet10_s_inline_simd_color_loopdiv_r3_proposed1.F
	$(FC) $(FCFLAGS) $^ -c -o $@

main.o: main.F
	$(FC) $(FCFLAGS) $^ -c -o $@

sub_data_io.o: sub_data_io.F
	$(FC) $(FCFLAGS) $^ -c -o $@

$(APP).x: main.o $(APP).o sub_data_io.o common/src/report.o
	$(FC) $(FCFLAGS) $^ -o $@


scop-detect:
	scop-detector . -e f -a flang $(FCFLAGS)

