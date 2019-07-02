#----------------------------------------------------------------------------
#   Makefile for MPI_PVODE.
#     $Id: Makefile,v 1.1.1.1 2008/09/04 05:54:22 dudson Exp $
#----------------------------------------------------------------------------

BOUT_TOP = ../../..

CINCLUDES = -I../include/pvode
FLAGS    = $(CXXFLAGS) $(CINCLUDES)

# PVODE library file
LIBF = ../lib/libpvode.a

# Source files
SRC.C = cvode.cpp nvector.cpp llnlmath.cpp cvspgmr.cpp spgmr.cpp iterativ.cpp cvdiag.cpp \
       smalldense.cpp

# Object files
ifdef NEED.IL
  OBJ.IL = $(SRC.C:%.cpp=obj/%.il)
endif

OBJS = $(SRC.C:%.cpp=obj/%.o) $(OBJ.IL)

all: $(LIBF)

$(LIBF): $(OBJS)
	$(AR) $(ARFLAGS) $(LIBF) $(OBJS)

obj:
	$(MKDIR) obj

obj/%.o: obj %.cpp
	$(CXX) $(FLAGS) -c -o $@ $*.cpp

clean::
	@rm -rf obj $(LIBF)
