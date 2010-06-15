#----------------------------------------------------------------------------
#   Makefile for MPI_PVODE preconditioners.
#     $Id: Makefile,v 1.1.1.1 2008/09/04 05:54:22 dudson Exp $
#----------------------------------------------------------------------------

BOUT_TOP = ../..
include $(BOUT_TOP)/make.config

CINCLUDES = -I$(PVODE)/include
FLAGS    = $(CFLAGS) $(CINCLUDES)

# PVODE library file
LIBF = ../lib/libpvpre.a

# Source files
SRC.C = pvbbdpre.c band.c

# Object files
ifdef NEED.IL
  OBJ.IL = $(SRC.C:%.c=obj/%.il)
endif

OBJS = $(SRC.C:%.c=obj/%.o) $(OBJ.IL)

all: $(LIBF)

$(LIBF): $(OBJS)
	$(AR) rcv $(LIBF) $(OBJS)
	$(RANLIB) $(LIBF)

obj:
	$(MKDIR) obj
		
obj/%.o: obj %.c %.h
	$(CC) $(FLAGS) -c -o $@ $*.c

clean::
	@rm -rf obj $(LIBF)
