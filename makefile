
BOUT_TOP  = .

DIRS      = src

TARGET   ?= libfast

include make.config

SOLIB=lib/libbout++.so.$(BOUT_VERSION)
STATLIB=lib/libbout++.a

shared-lib: $(SOLIB)

static-lib: $(STATLIB)


$(STATLIB): $(OBJ)
	@echo "  Creating static lib"
	$(AR) $(ARFLAGS) $@ $(OBJ)

$(SOLIB): $(OBJ)
	@echo "  Creating shared lib"
	@echo $(BOUT_FLAGS) | grep -i pic &>/dev/null || (echo "not compiled with PIC support - reconfigure with --enable-shared" ;exit 1)
	@$(RM) $(BOUT_TOP)/lib/*.so*
	@$(CXX) -shared -Wl,-soname,libpvode.so.1.0.0 -o $(BOUT_TOP)/lib/libpvode.so.1.0.0 -L $(BOUT_TOP)/lib -Wl,--whole-archive -lpvode -Wl,--no-whole-archive
	@cd lib ; ln -s libpvode.so.1.0.0 libpvode.so
	@$(CXX) -shared -Wl,-soname,libpvpre.so.1.0.0 -o $(BOUT_TOP)/lib/libpvpre.so.1.0.0 -L $(BOUT_TOP)/lib -Wl,--whole-archive -lpvpre -Wl,--no-whole-archive
	@cd lib ; ln -s libpvpre.so.1.0.0 libpvpre.so
	@$(CXX) -shared -Wl,-soname,libbout++.so.$(BOUT_VERSION) -o $(SOLIB) $(OBJ)
	@cd lib ; ln -s libbout++.so.$(BOUT_VERSION) libbout++.so

######################################################################
# Tests
######################################################################

check-unit-tests: libfast
	@export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH}; $(MAKE) --no-print-directory -C tests/unit check

check-mms-tests: libfast
	@cd tests/MMS; export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH} ; ./test_suite

check-integrated-tests: libfast
	@cd tests/integrated; LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH}./test_suite_make
	@cd tests/integrated; export LD_LIBRARY_PATH=${PWD}/lib:${LD_LIBRARY_PATH} ; \
		PYTHONPATH=${PWD}/tools/pylib/:${PYTHONPATH} ./test_suite


check: check-unit-tests check-integrated-tests check-mms-tests

build-check-unit-tests: libfast
	@$(MAKE) --no-print-directory -C tests/unit

build-check-mms-tests: libfast
	$(MAKE) --no-print-directory -C tests/MMS

build-check-integrated-tests: libfast
	$(MAKE) --no-print-directory -C tests/integrated


build-check: build-check-integrated-tests build-check-mms-tests build-check-unit-tests
