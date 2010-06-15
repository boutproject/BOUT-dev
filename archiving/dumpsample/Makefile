
CC = c++
LD = c++

CFLAGS = -Wall -g

PDB_PATH=/hwdisks/home/bd512/local/

INCLUDE = -I$(PDB_PATH)/include
LIBS = -lm -L$(PDB_PATH)/lib -lpdb -lpml -lscore 

TARGET = pdbsample
OBJ = pdbsample.o

.PHONY:all
all: $(TARGET)

$(TARGET): $(OBJ) Makefile
	$(LD) -o $(TARGET) $(OBJ) $(LIBS)

$(OBJ): %.o: %.cpp Makefile
	$(CC) $(CFLAGS) -c $< -o $@ $(INCLUDE)

.PHONY:clean
clean:
	rm -f $(OBJ) $(TARGET)

.PHONY:force
force: clean all
