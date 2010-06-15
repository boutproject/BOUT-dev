#----------------------------------------------------------------------------
#          Makefile for the EFITXS package
#
#----------------------------------------------------------------------------

CPUNAME = $(shell uname -n)

ifneq ($(findstring gre,$(CPUNAME)),)  #grendel.llnl.gov
  PDB_PATH = /usr/local/pact/pact04_05_11/
  IDL_DIR = /usr/local/rsi/idl
  CC = cc
  FF = pgf90
  LINKER = pgf90
endif

ifneq ($(findstring sma,$(CPUNAME)),)  #smaug.llnl.gov
  #use /afs/localcell/usr/rsi/idl70/bin/idl
  PDB_PATH = /afs/localcell/usr/pact/@sys/pact07_07_18/
  IDL_DIR = /afs/localcell/usr/rsi/idl70
  CC = cc
  FF = ifort ##f95
  LINKER = ifort ##f95
endif

ifneq ($(findstring hre,$(CPUNAME)),)  #hrethric.llnl.gov
  PDB_PATH = /usr/local/pact/pact_04_05_11/
  IDL_DIR = /usr/local/rsi/idl
  CC = cc
  FF = pgf90
  LINKER = pgf90
endif

ifneq ($(findstring jac,$(CPUNAME)),)  #jacquard.nersc.gov
  PDB_PATH =/usr/common/homes/u/umansky/OldPACT/PACT_jm/pact06_08_28/pact/dev/lnx-2.3-o/
  IDL_DIR = /usr/common/usg/idl/idl_6.3
  CC = cc
  FF = f95
  LINKER = f95
endif

ifneq ($(findstring sausage,$(CPUNAME)),)  # sausage at York
  PDB_PATH=/hwdisks/home/bd512/local/
  IDL_DIR=/hwdisks/sfw/idl/idl6.3/idl/
  CC=gcc
  FF=gfortran
  LINKER=gfortran
endif


CFILES = efitxs_wrap.c
FFILES = efitxs.f90

COBJECTS = efitxs_wrap.o
FOBJECTS = efitxs.o
EXECFIL = efitxs.so

X_LIBS =
LIBS =

CFLAGS = -I$(IDL_DIR)/external -fPIC -D UNDERSCORE
FFLAGS = -I$(IDL_DIR)/external -fPIC

LDFLAGS = -shared
X_LD_FLAGS = -I$(IDL_DIR)/external
EXPORT_CFLG = 
F_LD_POST= 



OBJECTS = $(COBJECTS) $(FOBJECTS)


$(EXECFIL) : $(OBJECTS)
	$(LINKER)  $(X_LD_FLAGS) $(LDFLAGS) $(EXPORT_FLAGS) -o $(EXECFIL) \
	$(OBJECTS) $(LIBS)
	$(F_LD_POST)



efitxs_wrap.o: efitxs_wrap.c
	$(CC) $(CFLAGS) -c efitxs_wrap.c

efitxs.o: efitxs.f90
	$(FF) $(FFLAGS) -c efitxs.f90


clean: 
	rm -f *.o *.so core.* *.mod *.export
