# configuration -- change this to install in different directories,
#                  or if "make" doesn't work.
# prefix = /usr/local

#########################
#       TARGET          #
#########################

TARGET= flex_dyn

SUF90=f90
SUF77=f

.SUFFIXES: .f90 .f .o .mod

GSL_prefix = /usr/local

#########################
#      COMPILATION      #
#########################

FAD	= mpif90
F90	= mpif90
F77	= mpif77


CC = gcc
FC = gfortran

FFLAGS  = -Wall -Wno-unused -Wno-tabs  #-r8 -O4 -openmp # -fpe3 -parallel  #-traceback #-ftrapuv -check uninit -traceback #  -g -fpe3 # -traceback -debug all

# -zero -fpe0  -CB  -O0  -g3 -debug extended -ftrapuv -check all #-parallel # -check
# -openmp #-check

CFLAGS = -fPIC -Wall -Wno-unused -g -I$(GSL_prefix)/include
LFLAGS =  #-L$(GSL_prefix)/lib -lgsl -lgslcblas -lm
LIBS = -ldl -lstdc++



#export:
#    ar rvs pcestimate.a *.o

SRCS = settings.o common_utils.o solver_utils.o dim_flexible_multi_body_dyn.o main.o

OBJS =  ${SRCS:.$(SUF)=.o}

all:  $(TARGET)

$(TARGET): $(OBJS) 
	$(F90) $(FFLAGS) -o $(TARGET) $(OBJS)  $(LFLAGS) -Wl,-rpath=.
	@echo " ----------- ${TARGET} created ----------- "

######################################
####### Compilation
######################################
%.o : %.mod

.$(SUF90).o:
	$(F90) $(FFLAGS) -c $<
.$(SUF77).o:
	$(F77) $(FFLAGS) -c $<

##################################
# Clean
##################################

clean:
	rm -f *.o *.so
