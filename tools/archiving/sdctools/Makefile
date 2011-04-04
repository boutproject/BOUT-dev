
SDC = ./sdclib/
PACT = /hwdisks/home/bd512/local/
IDL_DIR = /hwdisks/sfw/idl/idl6.3/idl/

CC = gcc
LD = gcc

CFLAGS = -Wall -O

INCLUDES = -I$(SDC) -I$(PACT)/include 
LIBS = $(SDC)/sdclib.o -L$(PACT)/lib -lpdb -lpml -lscore -lm


.PHONY: all

all: compress sdc2idl.so

compress: compress.o
	$(LD) -o compress compress.o $(LIBS) 

compress.o: compress.c
	$(CC) $(CFLAGS) -c compress.c $(INCLUDES)


sdc2idl.so: sdc2idl.o
	$(LD) -o sdc2idl.so sdc2idl.o $(SDC)/sdclib.so -shared -lm

sdc2idl.o: sdc2idl.c
	$(CC) -c $(CFLAGS) -fPIC -c sdc2idl.c -I$(SDC) -I$(IDL_DIR)/external


.PHONY: clean
clean:
	rm *.o compress sdc2idl.so
