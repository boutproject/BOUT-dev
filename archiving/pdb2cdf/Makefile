
CC = CC
LD = CC

CFLAGS = -g

PDB_PATH=$(HOME)/local/
CDF_PATH=$(HOME)/local/

INCLUDE = -I$(PDB_PATH)/include -I$(CDF_PATH)/include
LIBS = -lm -L$(PDB_PATH)/lib -lpdb -lpml -lscore -L$(CDF_PATH)/lib -lnetcdf_c++ -lnetcdf

TARGET = pdb2cdf
OBJ = pdb2cdf.o

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
