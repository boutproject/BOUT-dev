####################################################################
# Tests
####################################################################

BOUT_TOP = ..

BOUT_TEST_DIR = .

GTEST_DIR = $(BOUT_TOP)/googletest/googletest
GMOCK_DIR = $(BOUT_TOP)/googletest/googlemock

GTEST_SOURCES = $(GTEST_DIR)/src/gtest-all.cc

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CXXFLAGS += -isystem ${GTEST_DIR}/include -I${GTEST_DIR}

# Flags passed to the C++ compiler.
CXXFLAGS += -g -Wall -Wextra -pthread


# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

include $(BOUT_TOP)/make.config

# House-keeping build targets.

clean ::
	rm -f $(TESTS) gtest.a gtest_main.a *.o

# Builds gtest.a and gtest_main.a.

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
gtest-all.o : $(GTEST_SRCS_)
	@echo "  Compiling" $@
	@$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	@echo "  Compiling" $@
	@$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	@echo "  Linking" $@
	@$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	@echo "  Linking" $@
	@$(AR) $(ARFLAGS) $@ $^

$(TARGET): makefile $(BOUT_TOP)/make.config $(OBJ) $(SUB_LIBS) gtest_main.a
	@echo "  Linking" $(TARGET)
	@$(LD) $(LDFLAGS) -o $(TARGET) $(OBJ) $(BOUT_LIBS) $(SUB_LIBS) gtest_main.a

TEST_SOURCES = $(shell find $(BOUT_TEST_DIR) -type f -name "test_*.cxx" 2> /dev/null)
TEST_OBJECTS = $(TEST_SOURCES:%.cxx=%.o)

all_tests: makefile $(BOUT_TOP)/make.config $(OBJ) $(SUB_LIBS) $(TEST_OBJECTS) gtest_main.a
	@echo "  Linking tests"
	@$(LD) $(LDFLAGS) -o $@ $(OBJ) $(TEST_OBJECTS) gtest_main.a $(BOUT_LIBS) $(SUB_LIBS)

check: all_tests
	./all_tests

%.o: $(BOUT_TOP)/make.config %.cxx
	@echo "  Compiling " $(@:.o=.cxx)
	@$(CXX) $(BOUT_INCLUDE) $(BOUT_FLAGS) -c $(@:.o=.cxx) -o $(@)
