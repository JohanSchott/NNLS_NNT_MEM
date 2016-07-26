# Modules needed on matter:
# module load intel/12.1.4  mkl/10.3.10.319
# Modules needed on India:
# module load gcc/default  atlas/3.9.33/gcc

OBJS =  fz.o settingm.o matlab.o openm.o savem.o leastsquare.o tikhonovm.o asym.o eigenm.o maxentm.o input.o     smearm.o     readinpm.o acm.o  saveacm.o     main.o
PROG =  ./LSTMEM

#LIBRARIES on India:
LIBS = -L${LIBRARY_PATH}  -llapack -lf77blas -lcblas -latlas -L${ZGELSQUAD} -lzgelsquad -L${DGELSDQUAD} -ldgelsdquad
#LIBRARIES on matter:
#LIBS = -L$(ZGELSQUAD)  -lzgelsquad  -Nmkl -mkl=sequential
#LIBRARIES on local computer:
#LIBS =  -L${ZGELSQUAD} -lzgelsquad -L${DGELSDQUAD} -ldgelsdquad -L${LAPACKPATH} -llapack -lrefblas

# Compiler on matter:
#F90 = ifort 
#F90FLAGS =  -O2
# Compiler on India:
F90 = gfortran 
F90FLAGS =  -O2 # -fopenmp  -g -Wall -fbounds-check -Waliasing -fdump-core  -Wreturn-type 
# Compiler on local computer
#F90 = gfortran
#F90FLAGS = -O2 -g #-Wall -fbounds-check -Waliasing -fdump-core  -Wreturn-type

LDFLAGS = $(LIBS)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod


.SUFFIXES: $(SUFFIXES) .f90 .F90 .f .cpp

.f.o:
	$(F77) $(F77FLAGS) -c $<
.f90.o:
	$(F90) $(F90FLAGS) -c $<
.F90.o:
	$(F90) $(F90FLAGS) -c $<
.f90.mod:
	$(F90) $(F90FLAGS) -c $<
.F90.mod:
	$(F90) $(F90FLAGS) -c $<

.mod.o:
	$(F90) $(F90FLAGS) -c $<

.cpp.o:
	$(CC) $(INCLUDE) $(LIBS) -c $<

