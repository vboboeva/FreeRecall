################################################################################
#
# Makefile to compile and link C programs
#
# Version valid for Linux machines
#
# "make" compiles and links the specified main programs and modules
# using the specified libraries (if any), and produces the executables
# 
# "make clean" removes all files created by "make"
#
# Dependencies on included files are automatically taken care of
#
################################################################################

all: rmxeq mkdep mkxeq rmpmk
.PHONY: all

CC= gcc-7

# main programs and required modules 

MAIN = main #main_mpi
# test_discrete

RANDOM = ranlxs ranlxd gauss discrete

START = start utils

FUNCTIONS = functions #functions_mpi
#generate_patterns

MODULES = $(RANDOM) $(START) $(FUNCTIONS)


# search path for modules

MDIR = ../modules

VPATH = $(MDIR)/random:$(MDIR)/start:$(MDIR)/functions


# additional include directories

INCPATH = ../include  /usr/include
# /usr/include/gsl /usr/include/igraph  /usr/local/include/igraph  /u/shared/modules/libraries


# additional libraries to be included 

LIBS = m gomp
# gsl gslcblas igraph

LIBPATH = /usr/lib/gcc/x86_64-linux-gnu/7/


# scheduling and optimization options (such as -DSSE -DSSE2 -DP4)
 
CFLAGS =  -ftree-vectorize -mavx -ffast-math -O3 -std=gnu99 # -Wpedantic -Wall   -fopenmp -fopenmp-simd -fopt-info-vec-omp-optimized-missed=opt_log
#-lm -lgsl -lgslcblas
#-Wall 
#-lcblas
 

############################## do not change ###################################

SHELL=/bin/bash

PGMS= $(MAIN) $(MODULES)

INCDIRS = $(addprefix -I,$(INCPATH))

OBJECTS = $(addsuffix .o,$(MODULES))

LDFLAGS = $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(LIBS))

-include $(addsuffix .d,$(PGMS))


# rule to make dependencies

$(addsuffix .d,$(PGMS)): %.d: %.c Makefile
	@ $(CC) -MM $(INCDIRS) $< -o $@


# rule to compile source programs

$(addsuffix .o,$(PGMS)): %.o: %.c Makefile
	$(CC) $< -c $(CFLAGS) $(INCDIRS) -o $@


# rule to link object files

$(MAIN): %: %.o $(OBJECTS) Makefile
	$(CC) $< $(OBJECTS) $(CFLAGS) $(LDFLAGS) -o $@


# produce executables

mkxeq: $(MAIN)


# remove old executables and old error log file

rmxeq:
	@ -rm -f $(MAIN); \
        echo "delete old executables"		


# make dependencies

mkdep:  $(addsuffix .d,$(PGMS))
	@ echo "generate tables of dependencies"


# clean directory and executables only remain

rmpmk:
	@ -rm -rf *.d *.o .tmp


# clean directory 

clean:
	@ -rm -rf *.d *.o .tmp $(MAIN)
.PHONY: clean

################################################################################
