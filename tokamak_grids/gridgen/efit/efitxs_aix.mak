#----------------------------------------------------------------------------
#          Makefile for the EFITXS package
#
#----------------------------------------------------------------------------

#-for seaborg-
IDL_DIR = /usr/common/graphics/idl/idl

CC = cc
FF = f90
LINKER = cc

CFILES = efitxs_wrap.c
FFILES = efitxs.f90

COBJECTS = efitxs_wrap.o
FOBJECTS = efitxs.o
EXECFIL = efitxs.so

X_LIBS = -L$(IDL_DIR)/bin/bin.ibm -bI:$(IDL_DIR)/external/idl.export
LIBS = $(X_LIBS) -lm

CFLAGS = -I$(IDL_DIR)/external -q32 -D_THREAD_SAFE
FFLAGS = -I$(IDL_DIR)/external -q32

LDFLAGS = -bnoquiet
X_LD_FLAGS = -I$(IDL_DIR)/external -bM:SRE -bnoentry -q32 
EXPORT_CFLG = -bE:
F_LD_POST= -lxlf -lxlf90

#-this is needed for every function called from IDL
READ_DATA_G_EXPORT = read_data_g.export
SAVE_DATA_G_EXPORT = save_data_g.export
READ_DATA_A_EXPORT = read_data_a.export
SAVE_DATA_A_EXPORT = save_data_a.export
GET_DIMS_G_EXPORT = get_dims_g.export
GET_DATA_G_EXPORT = get_data_g.export
PUT_DATA_G_EXPORT = put_data_g.export
GET_DATA_A_EXPORT = get_data_a.export
PUT_DATA_A_EXPORT = put_data_a.export
WRITE_GRIDUE_EXPORT = write_gridue.export

#-list of all functions called by IDL
AIX_EXPORT = read_data_g.export save_data_g.export read_data_a.export save_data_a.export get_dims_g.export \
	get_data_g.export put_data_g.export get_data_a.export put_data_a.export write_gridue.export


OBJECTS = $(COBJECTS) $(FOBJECTS)


$(EXECFIL) : $(OBJECTS) $(AIX_EXPORT)
	$(LINKER)  $(F_LD_FLAGS) $(X_LD_FLAGS) $(EXPORT_FLAGS) -o $(EXECFIL) \
	$(OBJECTS) $(LIBS)\
	$(EXPORT_CFLG)$(READ_DATA_G_EXPORT) \
	$(EXPORT_CFLG)$(SAVE_DATA_G_EXPORT) \
	$(EXPORT_CFLG)$(READ_DATA_A_EXPORT) \
	$(EXPORT_CFLG)$(SAVE_DATA_A_EXPORT) \
	$(EXPORT_CFLG)$(GET_DIMS_G_EXPORT) \
	$(EXPORT_CFLG)$(GET_DATA_G_EXPORT) \
	$(EXPORT_CFLG)$(PUT_DATA_G_EXPORT) \
	$(EXPORT_CFLG)$(GET_DATA_A_EXPORT) \
	$(EXPORT_CFLG)$(PUT_DATA_A_EXPORT) \
	$(EXPORT_CFLG)$(WRITE_GRIDUE_EXPORT) \
	$(F_LD_POST)



efitxs_wrap.o: efitxs_wrap.c
	$(CC) $(CFLAGS) -c efitxs_wrap.c

efitxs.o: efitxs.f90
	$(FF) $(FFLAGS) -c efitxs.f90


#-needed for every function called by IDL
read_data_g.export:
	echo "read_data_g" > read_data_g.export 

save_data_g.export:
	echo "save_data_g" > save_data_g.export 

read_data_a.export:
	echo "read_data_a" > read_data_a.export 

save_data_a.export:
	echo "save_data_a" > save_data_a.export 

get_dims_g.export:
	echo "get_dims_g" > get_dims_g.export 

get_data_g.export:
	echo "get_data_g" > get_data_g.export 

put_data_g.export:
	echo "put_data_g" > put_data_g.export 

get_data_a.export:
	echo "get_data_a" > get_data_a.export 

put_data_a.export:
	echo "put_data_a" > put_data_a.export 

write_gridue.export:
	echo "write_gridue" > write_gridue.export

clean: 
	rm -f *.o *.so core.* *.mod *.export
