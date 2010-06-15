TARGET = 2fluid
SRC = 2fluid.cpp

#############################################
# The rest of this is the same for all tests

ifndef BOUT_TOP
	BOUT_TOP = ../..
endif

include $(BOUT_TOP)/make.config

OBJ = $(SRC:%.cpp=%.o)

$(TARGET): $(BOUT_TOP)/lib/libbout++.a $(OBJ) $(BOUT_TOP)/make.config Makefile
	$(LD) -o $(TARGET) $(OBJ) $(BOUT_LIBS)

%.o: %.cpp $(BOUT_TOP)/make.config Makefile
	$(CXX) $(CFLAGS) -c $< -o $@ $(BOUT_INCLUDE)

.PHONY: clean
clean:
	rm $(OBJ) $(TARGET)
