
CFLAGS = `bout-config --cflags`
LIBS = `bout-config --libs`
LD = `bout-config --ld`

test_cbout: test_cbout.o
	$(LD) test_cbout.o $(CFLAGS) $(LIBS)

test_cbout.o: test_cbout.c
	gcc -c test_cbout.c -I../../../include/

