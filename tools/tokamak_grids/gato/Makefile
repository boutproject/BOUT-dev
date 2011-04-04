
CC = gcc
LD = gcc

CFLAGS = -Wall

PDB_PATH=/usr/local/pact/pact

INCLUDE = -I$(PDB_PATH)/include
LIBS = -lm -L$(PDB_PATH)/lib -lpdb -lpml -lscore

TARGET = gato2pdb
OBJ = gato2pdb.o read_gato.o

.PHONY:all
all: $(TARGET)

$(TARGET): $(OBJ)
	$(LD) -o $(TARGET) $(OBJ) $(LIBS)

$(OBJ): %.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@ $(INCLUDE)
