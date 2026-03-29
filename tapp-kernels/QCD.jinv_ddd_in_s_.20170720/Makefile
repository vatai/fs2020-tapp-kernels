CC=clang
CXX=clang++
CFLAGS = -fopenmp -Icommon/include -DRDC -DVLENS=16 -DEML_LIB -DPREFETCH -D_CHECK_SIM -DTARGET_JINV
CXXFLAGS = -std=gnu++11  $(CFLAGS)
LDLIBS= -lrand

prog.x: bicgstab_precdd_s.o clover_s.o ddd_in_s.o ddd_out_s.o main.o qws.o static_solver.o common/src/report.o
	$(CXX) $^ -fopenmp -o $@
