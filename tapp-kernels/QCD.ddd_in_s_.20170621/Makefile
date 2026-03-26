CC=clang
CXX=clang++
CFLAGS = -fopenmp -Icommon/include -DRDC -DVLENS=16 -DEML_LIB -DPREFETCH -D_CHECK_SIM -DTARGET_JINV
CXXFLAGS = -std=gnu++11  $(CFLAGS)
LDLIBS= -lrand

prog.x: bicgstab_dd_mix.o bicgstab_precdd_s.o clover_s.o ddd_in_s.o main.o qws.o common/src/report.o
	$(CXX) $^ -fopenmp -o $@

scop-detect:
	scop-detector . -e c -a clang++ -std=gnu++11  -DRDC -DV512 -Icommon/include -fopenmp
