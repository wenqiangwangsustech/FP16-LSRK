#================================================================
#   Copyright (C) 2021 Sangfor Ltd. All rights reserved.
#   
#   File Name：Makefile
#   Author: Wenqiang Wang
#   Created Time:2021-10-30
#   Discription:
#
#================================================================

GPU_CUDA := ON


LayerMedium := #ON


FREE_SURFACE := ON
PML := ON
SOLVE_DISPLACEMENT := ON

OldEq := ON

ifdef OldEq
FLOAT16 := 
#old equation
SRCDIR := ./src/src_

#willimson Runge-Kutta
#SRCDIR := ./src/src_willimson
else
FLOAT16 := ON

#All physical variables
#SRCDIR := ./src/src_wave_con_med

#willimson Runge-Kutta and All physical variables
SRCDIR := ./src/src_willimson_wave_con

#contravariant jacobian and medium
#SRCDIR := ./src/src_con_med_fp


#willimson Runge-Kutta and contravariant jacobian and medium
#SRCDIR := ./src/src_willimson_con
endif






CCHOME := /usr
CUDAHOME := /public/software/cuda-10.0
MPIHOME := /public/software/openmpi-4.1.1-cuda.10
PROJHOME := /public/software/proj-8.1.0



CC := $(CCHOME)/bin/gcc -pipe

#General Compiler
ifdef GPU_CUDA
GC := $(CUDAHOME)/bin/nvcc -rdc=true -maxrregcount=127 -arch=sm_70 #-Xptxas=-v 
else
GC := $(CCHOME)/bin/g++ -pipe
endif


LIBS := -L$(CUDAHOME)/lib64 -lcudart -lcublas
INCS := -I$(CUDAHOME)/include 

LIBS += -L$(MPIHOME)/lib -lmpi
INCS += -I$(MPIHOME)/include 


LIBS += -L$(PROJHOME)/lib -lproj
INCS += -I$(PROJHOME)/include  



OBJDIR := ./obj
BINDIR := ./bin


CFLAGS := -c -O2 
LFLAGS := -O2

GCFLAGS := 

ifdef GPU_CUDA
#LFLAGS += -Xptxas=-v 

#LFLAGS += -arch=sm_70 -rdc=true -Xptxas=-v 
#GCFLAGS += --fmad=false 
GCFLAGS += -x cu
endif

vpath

vpath % $(SRCDIR)
vpath % $(OBJDIR)
vpath % $(BINDIR)


DFLAGS_LIST := GPU_CUDA FLOAT16 \
			   FREE_SURFACE PML SOLVE_DISPLACEMENT LayerMedium


DFLAGS := $(foreach flag,$(DFLAGS_LIST),$(if $($(flag)),-D$(flag)))


OBJS := cjson.o printInfo.o create_dir.o readParams.o \
		init_gpu.o cpu_Malloc.o init_grid.o init_MPI.o \
		run.o \
		coord.o terrain.o medium.o dealMedium.o crustMedium.o calc_CFL.o\
		contravariant.o \
		wave_deriv.o wave_rk.o freeSurface.o \
		init_pml_para.o pml_deriv.o pml_freeSurface.o \
		propagate.o \
		pml_rk.o \
		singleSource.o multiSource.o\
		data_io.o \
		MPI_send_recv.o MPI_send_recv_jac.o \
		PGV.o station.o\
		main.o #MPI_send_recv_fp32.o




OBJS := $(addprefix $(OBJDIR)/,$(OBJS))


$(BINDIR)/CGFDM3D: $(OBJS)
	$(GC) $(LFLAGS) $(LIBS) $^ -o $@


$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(GC) $(CFLAGS) $(DFLAGS) $(GCFLAGS) $(INCS)  $^ -o $@


$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) $^ -o $@

clean:
	-rm $(OBJDIR)/* -rf
	-rm $(BINDIR)/* -rf
	-rm output -rf
