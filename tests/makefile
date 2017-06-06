####################################################################
# Tests
####################################################################

BOUT_TOP = ..

BOUT_TEST_DIR = .

TESTS = gtest

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

TESTS = msg_stack_test

include $(BOUT_TOP)/make.config

all : $(TESTS)

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
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

# Builds a sample test.  A test should link with either gtest.a or
# gtest_main.a, depending on whether it defines its own main()
# function.
# sample1.o : $(USER_DIR)/sample1.cc $(USER_DIR)/sample1.h $(GTEST_HEADERS)
# 	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/sample1.cc

# sample1_unittest.o : $(USER_DIR)/sample1_unittest.cc \
#                      $(USER_DIR)/sample1.h $(GTEST_HEADERS)
# 	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(USER_DIR)/sample1_unittest.cc

# sample1_unittest : sample1.o sample1_unittest.o gtest_main.a
# 	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -lpthread $^ -o $@

SOURCEC = msg_stack_test.cxx

$(TARGET): makefile $(BOUT_TOP)/make.config $(OBJ) $(SUB_LIBS) gtest_main.a
	@echo "  Linking" $(TARGET)
	@$(LD) $(LDFLAGS) -o $(TARGET) $(OBJ) $(BOUT_LIBS) $(SUB_LIBS) gtest_main.a

