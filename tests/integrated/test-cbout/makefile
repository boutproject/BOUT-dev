
CFLAGS = $(shell bout-config --cflags)
LIBS = $(shell bout-config --libs)
LD=$(shell bout-config --ld)

test_cbout: test_cbout.o
	$(LD) -o $@ $^ $(CFLAGS) $(LIBS)

test_cbout.o: test_cbout.c 
	gcc -c test_cbout.c -I../../../include/

