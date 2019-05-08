# Example makefile. Sources are listed in SOURCEC
# and optionally a TARGET executable name can be given
# It requires that bout-config is in $PATH

SOURCEC = test.cxx

############

OBJ = $(SOURCEC:%.cxx=%.o)
TARGET ?= $(SOURCEC:%.cxx=%)

# Use the bout-config script to get compiler, flags
CXX:=$(shell bout-config --cxx)

ifneq "$(.SHELLSTATUS)" "0"
which=$(shell which bout-config 2>&1 )
$(error Is bout-config in $$PATH? $(which))
endif

CFLAGS:=$(shell bout-config --cflags)
LD:=$(shell bout-config --ld)
LDFLAGS:=$(shell bout-config --libs)

$(TARGET): makefile $(OBJ)
	@echo "  Linking" $(TARGET)
	@$(LD) -o $(TARGET) $(OBJ) $(LDFLAGS)

%.o: %.cxx
	@echo "  Compiling " $(@F:.o=.cxx)
	@$(CXX) $(CFLAGS) -c $(@F:.o=.cxx) -o $@

.PHONY: clean
clean:
	rm $(OBJ) $(TARGET)

