BOUT_TOP	= ../../..
include $(BOUT_TOP)/make.config

TARGET = pdb2cdf
OBJ = pdb2cdf.o

.PHONY:all
all: $(TARGET)

$(TARGET): $(OBJ) Makefile
	$(LD) -o $(TARGET) $(OBJ) $(EXTRA_LIBS)

$(OBJ): %.o: %.cxx Makefile
	$(CC) $(CFLAGS) -c $< -o $@ $(EXTRA_INCS)
