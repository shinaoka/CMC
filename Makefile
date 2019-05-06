#For Xeon with SSE2
#CC=mpicc
#CPC=mpiCC
#CC=openmpic
#CPC=openmpic++
#CC=icc
#CPC=icpc
#FFLAGS = -O3 -DDSFMT_MEXP=521 -DHAVE_SSE2 -D_mpi_use -D_utime
#FFLAGS = -O3 -DDSFMT_MEXP=521 -DHAVE_SSE2 -D_mpi_use -xSSE4.1 -axSSE4.1 -D_utime
#FFLAGS = -O3 -DDSFMT_MEXP=521 -D_mpi_use -D_utime
#FFLAGS = -g -DDSFMT_MEXP=521 -DHAVE_SSE2 -D_utime

#Brutus
#CC=mpicc
#CPC=mpic++
#FFLAGS = -O3 -DDSFMT_MEXP=521 -DHAVE_SSE2 -D_mpi_use -D_utime

#For System A
#CC=mpixcc
#FFLAGS = -Os -64 -noparallel  -DDSFMT_MEXP=521 -lm -D_mpi_use -D_ising

#For System B
#CC=icc
#CPC=icpc
#FFLAGS = -O3 -ipo -DDSFMT_MEXP=521 -lmpi -D_mpi_use 
#FFLAGS = -O3 -ipo -DDSFMT_MEXP=521 

#gcc
#CC=gcc
#CPC=g++
#FFLAGS = -O3 -DDSFMT_MEXP=521 -DHAVE_SSE2 -msse2 

#gcc4.4+openmpi (MacPorts)
#CC=/opt/local/bin/openmpicc
#CPC=/opt/local/bin/openmpic++
#FFLAGS = -O3 -DDSFMT_MEXP=521 -DHAVE_SSE2 -D_mpi_use -D_utime -msse4.2
#FFLAGS = -g -DDSFMT_MEXP=521 -DHAVE_SSE2 -D_utime

#Monch
CC=mpicc
CPC=mpicxx
FFLAGS = -O3 -DDSFMT_MEXP=521 -DHAVE_SSE2 -D_mpi_use -D_utime -I/home/tool/opt/boost_1_60_0/include

#CC=mpicc
#CPC=icpc
#FFLAGS = -O3 -DDSFMT_MEXP=521 -DHAVE_SSE2 -D_utime -I$(HOME)/src/boost_1_55_0

LFLAGS = $(FFLAGS)
LIB = 

.SUFFIXES :
.SUFFIXES : .o .c
.SUFFIXES : .o .cpp

.c.o:
	$(CC) -c $(FFLAGS) $< $(LIB)
.cpp.o:
	$(CPC) -c $(FFLAGS) $< $(LIB)

all: main.o mc.o util.o dSFMT.o phys_value.o loop_move.o ham.o phys_value_recorder.o mpi_util.o update.o
	$(CPC) $(FFLAGS) -o mcmain main.o mc.o util.o dSFMT.o phys_value.o loop_move.o ham.o phys_value_recorder.o mpi_util.o update.o $(LIB)

clean :
	rm -f *.o mcmain *.core

main.o: main.cpp mc.o util.o dSFMT.o mpi_util.o
mc.o: mc.cpp util.o phys_value.cpp exchange.cpp update.cpp loop_move.o ham.o dSFMT.o phys_value_recorder.o update.o
phys_value.o : phys_value.cpp
loop_move.o : loop_move.cpp dSFMT.o
ham.o : ham.cpp
dSFMT.o : dSFMT.cpp
update.o : update.cpp util.o
