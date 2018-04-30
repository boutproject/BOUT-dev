
BOUT_TOP  = .

DIRS      = src

TARGET   ?= libfast

include make.config

shared: libfast
	@echo "Creating libbout++.so"
	@echo $(BOUT_FLAGS) | grep -i pic > /dev/null 2>&1 || (echo "not compiled with PIC support - reconfigure with --enable-shared" ;exit 1)
	@#$(CXX) -shared -o $(LIB_SO) $(shell find $(BOUT_TOP)/src -name \*.o -type f -print 2> /dev/null) -L $(BOUT_TOP)/lib -Wl,--whole-archive -lpvode -lpvpre  -Wl,--no-whole-archive
	@$(RM) $(BOUT_TOP)/lib/*.so*
	@$(CXX) -shared -Wl,-soname,libbout++.so.$(BOUT_VERSION) -o $(LIB_SO).$(BOUT_VERSION) $(shell find $(BOUT_TOP)/src -name \*.o -type f -print 2> /dev/null)
	@$(CXX) -shared -Wl,-soname,libpvode.so.1.0.0 -o $(BOUT_TOP)/lib/libpvode_.so -L $(BOUT_TOP)/lib -Wl,--whole-archive -lpvode -Wl,--no-whole-archive
	@$(CXX) -shared -Wl,-soname,libpvpre.so.1.0.0 -o $(BOUT_TOP)/lib/libpvpre_.so -L $(BOUT_TOP)/lib -Wl,--whole-archive -lpvpre -Wl,--no-whole-archive
	@mv $(BOUT_TOP)/lib/libpvode_.so $(BOUT_TOP)/lib/libpvode.so.1.0.0
	@mv $(BOUT_TOP)/lib/libpvpre_.so $(BOUT_TOP)/lib/libpvpre.so.1.0.0
	@ln -s libbout++.so.$(BOUT_VERSION) $(LIB_SO)
	@ln -s libpvode.so.1.0.0 lib/libpvode.so
	@ln -s libpvpre.so.1.0.0 lib/libpvpre.so


######################################################################
# Tests
######################################################################

check-unit-tests: libfast
	@export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH}; $(MAKE) --no-print-directory -C tests/unit check

check-mms-tests: libfast
	@cd tests/MMS; export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH} ; \
		PYTHONPATH=${PWD}/tools/pylib/:${PYTHONPATH} ./test_suite_make
	@cd tests/MMS; export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH} ; \
		PYTHONPATH=${PWD}/tools/pylib/:${PYTHONPATH} ./test_suite

check-mms-tests-all: libfast
	@cd tests/MMS; export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH} ; \
		PYTHONPATH=${PWD}/tools/pylib/:${PYTHONPATH} ./test_suite_make --set-bool slow_tests=True
	@cd tests/MMS; export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH} ; \
		PYTHONPATH=${PWD}/tools/pylib/:${PYTHONPATH} ./test_suite --set-bool slow_tests=True

check-integrated-tests: libfast
	@cd tests/integrated; export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH} ; \
		PYTHONPATH=${PWD}/tools/pylib/:${PYTHONPATH} ./test_suite_make
	@cd tests/integrated; export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH} ; \
		PYTHONPATH=${PWD}/tools/pylib/:${PYTHONPATH} ./test_suite

check-integrated-tests-all: libfast
	@cd tests/integrated; export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH} ; \
		PYTHONPATH=${PWD}/tools/pylib/:${PYTHONPATH} ./test_suite_make --set-bool slow_tests=True
	@cd tests/integrated; export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH} ; \
		PYTHONPATH=${PWD}/tools/pylib/:${PYTHONPATH} ./test_suite --set-bool slow_tests=True


check: check-unit-tests check-integrated-tests check-mms-tests
check-all: check-unit-tests check-integrated-tests-all check-mms-tests-all

build-check-unit-tests: libfast
	@$(MAKE) --no-print-directory -C tests/unit

build-check-mms-tests: libfast
	$(MAKE) --no-print-directory -C tests/MMS

build-check-integrated-tests: libfast
	$(MAKE) --no-print-directory -C tests/integrated


build-check: build-check-integrated-tests build-check-mms-tests build-check-unit-tests
