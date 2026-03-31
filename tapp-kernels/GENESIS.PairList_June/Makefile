#for Intel
#FC	= gfortran
#FCFLAGS = -O3 -g -Wall  -ffree-line-length-none  -fallow-argument-mismatch -fopenmp
#LDFLAGS = -fopenmp
# for Fugaku
FC	= flang
#FCFLAGS = -O0 -g -Wall  -ffree-line-length-none -fbacktrace -fallow-argument-mismatch -fopenmp
FCFLAGS = -fopenmp
LDFLAGS = -fopenmp

OBJS= sub_data_io.o initialize_data.o kernel_pairlist_june.o  main.o
#OBJS_conv= sub_data_io.o convert_data.o

kernel : $(OBJS)
	$(FC) -o kernel sub_data_io.o initialize_data.o kernel_pairlist_june.o  main.o $(LDFLAGS)
#conv : $(OBJS_conv)
#	$(FC) -o conv sub_data_io.o convert_data.o $(LDFLAGS)
.SUFFIXES : .o .f90
.f90.o: 
	$(FC) $(FCFLAGS) -c $<
.PHONY: clean
clean:
	rm -f *.o *.mod kernel conv *.bak
