
CC = gcc
LD = gcc

CFLAGS = -Wall

PDB_PATH=/hwdisks/home/bd512/local/

INCLUDE = -I$(PDB_PATH)/include
LIBS = -lm -L$(PDB_PATH)/lib -lpdb -lpml -lscore

TARGET = elite2pdb
OBJ = elite2pdb.o

.PHONY:all
all: $(TARGET)

$(TARGET): $(OBJ)
	$(LD) -o $(TARGET) $(OBJ) $(LIBS)

$(OBJ): %.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@ $(INCLUDE)
