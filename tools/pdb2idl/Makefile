
include make.config

TARGET = pdb2idl.so

SRC = pdb2idl.c utils.c

OBJ = $(SRC:%.c=%.o)

.PHONY: all
all: $(TARGET)

$(TARGET): $(OBJ) make.config Makefile
	$(LD) -o $(TARGET) $(OBJ) $(LDFLAGS) $(LIB)

#### HEADER DEPENDENCIES ####
$(OBJ): utils.h

%.o: %.c make.config Makefile
	$(CC) $(CFLAGS) -c $< -o $@ $(INCLUDE)

.PHONY: clean
clean:
	@rm -f $(OBJ) $(TARGET)
