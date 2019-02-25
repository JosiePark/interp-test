CC = gfortran
LIBS = -L/usr/lib  -L/usr/local/lib -L/usr/lib/lapack -L/usr/lib/libblas

LINKS = -lnetcdf -lnetcdff -lfftw3  -llapack  -lblas
HOME_DIR = /home/josiepark/Project/PhD/CODE/INTERP/TESTS
#HOME_DIR = /home/josiepark/hpc-cluster/JOSIE_CABARET_test

INCLUDES = -I/usr/include  -I/usr/local/include -I$(build_DIR) -I$(HOME_DIR)/2WAVE

#CFLAGS = -O2 -g -heap-arrays -traceback -check all -fp-stack-check 
#NETCDF_VERSION=4.3.3
#CC = ifort

#LIBS = -L${MKL_HOME}/interfaces/fftw3xf -L/apps/netcdf/$(NETCDF_VERSION)/lib -L${MKL_HOME}/interfaces
#LINKS = -mkl -lfftw3xf_intel -lnetcdff -lnetcdf -llapack  -lblas

#HOME_DIR = /export111/work/jp1115/JOSIE_CABARET
NETCDF_DIR = $(HOME_DIR)/MODULES/NETCDF
LAGR_DIR = $(HOME_DIR)/MODULES/LAGR
#VAR_DIR = $(HOME_DIR)/src/VAR
#SRC_DIR = $(HOME_DIR)/src/TRANSPORT
#CAB_DIR = $(HOME_DIR)/src/CABARET
#SOLVER_DIR = $(HOME_DIR)/src/SOLVER
build_DIR = $(HOME_DIR)/BUILD
TWOWAVE_DIR = $(HOME_DIR)/2WAVE
JET_DIR = $(HOME_DIR)/JET
STOMMEL_DIR = $(HOME_DIR)/STOMMEL

#INCLUDES = -I${MKL_HOME}/include/fftw -I${MKL_HOME}/include -I/apps/netcdf/$(NETCDF_VERSION)/include -I$(build_DIR) -I$(SRC_DIR) 


