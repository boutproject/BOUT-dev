CC = gcc
LD = gcc

CFLAGS = -Wall

.PHONY: all

all: sdc_test sdclib.so

sdc_test: sdclib.o sdc_test.c
	$(CC) $(CFLAGS) -c sdc_test.c
	$(LD) -o sdc_test sdc_test.o sdclib.o -lm

sdclib.o: sdclib.c sdclib.h
	$(CC) $(CFLAGS) -c sdclib.c

sdclib.so: sdclib.c sdclib.h
	$(CC) $(CFLAGS) -fPIC -o sdclib.so -c sdclib.c

.PHONY: clean

clean:
	rm sdc_test sdc_test.o sdclib.o sdclib.so
