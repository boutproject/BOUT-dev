# Example makefile. Sources are listed in SOURCEC
# and optionally a TARGET executable name can be given

SOURCEC = test.cxx

############

OBJ = $(SOURCEC:%.cxx=%.o)
TARGET ?= $(SOURCEC:%.cxx=%)

# Use the bout-config script to get compiler, flags
CXX=`bout-config --cxx`
CFLAGS=`bout-config --cflags`
LD=`bout-config --ld`
LDFLAGS=`bout-config --libs`

$(TARGET): makefile $(OBJ)
	@echo "  Linking" $(TARGET)
	@$(LD) -o $(TARGET) $(OBJ) $(LDFLAGS)

%.o: %.cxx
	@echo "  Compiling " $(@F:.o=.cxx)
	@$(CXX) $(CFLAGS) -c $(@F:.o=.cxx) -o $@

.PHONY: clean
clean:
	rm $(OBJ) $(TARGET)

