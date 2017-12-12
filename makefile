
BOUT_TOP	= .

DIRS			= src

ifndef TARGET
TARGET=libfast
endif
# Add this to DIRS to have examples compiled
#examples

include make.config

shared: libfast
	@echo "Creating libbout++.so"
	@echo $(BOUT_FLAGS) | grep -i pic &>/dev/null || (echo "not compiled with PIC support - reconfigure with --enable-shared" ;exit 1)
	@#$(CXX) -shared -o $(LIB_SO) $(shell find $(BOUT_TOP)/src -name \*.o -type f -print 2> /dev/null) -L $(BOUT_TOP)/lib -Wl,--whole-archive -lpvode -lpvpre  -Wl,--no-whole-archive
	@$(CXX) -shared -o $(LIB_SO) $(shell find $(BOUT_TOP)/src -name \*.o -type f -print 2> /dev/null)
	@$(RM) $(BOUT_TOP)/lib/libpvode.so
	@$(RM) $(BOUT_TOP)/lib/libpvpre.so
	@$(CXX) -shared -o $(BOUT_TOP)/lib/libpvode_.so -L $(BOUT_TOP)/lib -Wl,--whole-archive -lpvode -Wl,--no-whole-archive
	@$(CXX) -shared -o $(BOUT_TOP)/lib/libpvpre_.so -L $(BOUT_TOP)/lib -Wl,--whole-archive -lpvpre -Wl,--no-whole-archive
	@mv $(BOUT_TOP)/lib/libpvode_.so $(BOUT_TOP)/lib/libpvode.so
	@mv $(BOUT_TOP)/lib/libpvpre_.so $(BOUT_TOP)/lib/libpvpre.so


######################################################################
# Tests
######################################################################

check-unit-tests:
	@export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH}; $(MAKE) --no-print-directory -C tests/unit check

check-mms-tests:
	@cd tests/MMS; export LD_LIBRARY_PATH=${PWD}/../../lib:${LD_LIBRARY_PATH} ; ./test_suite

check-integrated-tests:
	@cd tests/integrated; LD_LIBRARY_PATH=${PWD}/../../lib:${LD_LIBRARY_PATH}./test_suite_make
	@cd tests/integrated; export LD_LIBRARY_PATH=${PWD}/../../lib:${LD_LIBRARY_PATH} ; \
		PYTHONPATH=${PWD}/../../tools/pylib/:${PYTHONPATH} ./test_suite


check: check-unit-tests check-integrated-tests check-mms-tests
