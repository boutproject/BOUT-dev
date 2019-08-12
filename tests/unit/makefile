####################################################################
# Tests
####################################################################

BOUT_TOP = ../..

BOUT_TEST_DIR = .

GTEST_BASE = $(BOUT_TOP)/externalpackages/googletest
GTEST_DIR = $(GTEST_BASE)/googletest
GMOCK_DIR = $(GTEST_BASE)/googlemock

GTEST_SOURCES = $(GTEST_DIR)/src/gtest-all.cc

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CXXFLAGS += -isystem ${GTEST_DIR}/include -I${GTEST_DIR} -I${BOUT_TEST_DIR}
CXXFLAGS += -isystem ${GMOCK_DIR}/include -I${GMOCK_DIR}

# Flags passed to the C++ compiler.
CXXFLAGS += -g -Wall -Wextra -pthread


# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

all:  serial_tests
include $(BOUT_TOP)/make.config

LDFLAGS += -pthread

# Existence of this file indicates that GTEST has been downloaded
GTEST_SENTINEL = $(GTEST_BASE)/README.md

# House-keeping build targets.

clean ::
	@$(RM) -f gtest.a bout_test_main.a *.o $(TEST_OBJECTS)

# Builds gtest.a and bout_test_main.a.

# Download Google Test
$(GTEST_SENTINEL) :
	@echo "Downloading Google Test"
	git submodule update --init --recursive

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
#
# Note: If the GoogleTest submodule has not been checked out when
# this makefile is run, then the list of .cc and .h files is not
# correctly captured. Since we will probably not be modifying googletest
# this rule just depends on the sentinel file.
gtest-all.o : $(GTEST_SENTINEL)
	@echo "  Compiling" $@
	@$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(BOUT_FLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc

bout_test_main.o : $(GTEST_SENTINEL) $(BOUT_TEST_DIR)/bout_test_main.cxx
	@echo "  Compiling" $@
	@$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(BOUT_INCLUDE) $(BOUT_FLAGS) -c $(@:.o=.cxx)

gtest.a : gtest-all.o
	@echo "  Linking" $@
	@$(AR) $(ARFLAGS) $@ $^

bout_test_main.a : gtest-all.o bout_test_main.o
	@echo "  Linking" $@
	@$(AR) $(ARFLAGS) $@ $^

$(TARGET): makefile $(BOUT_TOP)/make.config $(OBJ) $(SUB_LIBS) bout_test_main.a
	@echo "  Linking" $(TARGET)
	@$(LD) $(LDFLAGS) -o $(TARGET) $(OBJ) $(BOUT_LIBS) $(SUB_LIBS) bout_test_main.a

TEST_SOURCES = $(shell find $(BOUT_TEST_DIR) -type f -name "test_*.cxx" 2> /dev/null)
TEST_OBJECTS = $(TEST_SOURCES:%.cxx=%.o)

$(TEST_OBJECTS): | $(GTEST_SENTINEL)

serial_tests: makefile $(BOUT_TOP)/make.config $(OBJ) $(LIB) checklib \
              $(SUB_LIBS) $(TEST_OBJECTS) bout_test_main.a
	@echo "  Linking tests"
	@$(LD) -o $@ $(TEST_OBJECTS) bout_test_main.a $(LDFLAGS) $(BOUT_LIBS) $(SUB_LIBS)

# Override this on the command line to print full output from the unit tests
QUIET ?= ./quiet.sh

# Note: This depends on GTEST_SENTINEL, which checks out Google Test if not present
#       The correct behaviour relies on the dependencies being checked in left-to-right order
#       The GTEST_SENTINEL dependency cannot be added to the %.o target, since this results in the %.o
#       target in make.config being used instead.
check: $(GTEST_SENTINEL) serial_tests
	@echo "Running unit test"
	@$(QUIET) ./serial_tests $(GTEST_ARGS)

# This is the same target as a make rule in make.config. Adding unmet dependencies here
# results in makefile reverting to the make.config version.
%.o: $(BOUT_TOP)/make.config %.cxx
	@echo "  Compiling " $(@:.o=.cxx)
	@$(CXX) $(BOUT_INCLUDE) $(BOUT_FLAGS) -c $(@:.o=.cxx) -o $(@)
