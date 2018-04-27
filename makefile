
BOUT_TOP  = .

DIRS      = src

TARGET   ?= libfast

all: main-target

include make.config

######################################################################
# Library
######################################################################

SOLIB=lib/libbout++.so.$(BOUT_VERSION)
STATLIB=lib/libbout++.a

shared-lib: $(SOLIB)

static-lib: $(STATLIB)

$(STATLIB): $(OBJ)
	@echo "  Creating static lib"
	@$(AR) $(ARFLAGS) $@ $(OBJ)

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

distclean:: clean clean-tests
# Removing the externalpackage installation. When we have more packages, need a better way
	@$(RM) -rf $(BOUT_TOP)/include/pvode
	@echo lib cleaned
	@$(RM) -rf $(BOUT_TOP)/lib/*
	-@$(RM) $(BOUT_TOP)/externalpackages/PVODE/lib/*.a
	-@$(RM) $(BOUT_TOP)/externalpackages/PVODE/source/obj/*.o
	-@$(RM) $(BOUT_TOP)/externalpackages/PVODE/precon/obj/*.o
	-@$(RM) -rf $(BOUT_TOP)/autom4te.cache make.config.{old,new}
	@echo externalpackages cleaned
	@echo autom4te.cache cleaned

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

clean-tests: clean-unit-tests clean-integrated-tests clean-mms-tests

clean-unit-tests:
	@echo "   tests/unit cleaned"
	@$(MAKE) --no-print-directory -C tests/unit clean

clean-integrated-tests:
	@echo "   tests/integrated cleaned"
	@$(MAKE) --no-print-directory -C tests/integrated clean

clean-mms-tests:
	@echo "   tests/MMS cleaned"
	@$(MAKE) --no-print-directory -C tests/MMS clean

####################################################################
# Documentation
####################################################################

MANUAL_DIR=$(BOUT_TOP)/manual

doxygen:
	@$(MAKE) -C $(MANUAL_DIR) doxygen

breathe-autogen:
	@$(MAKE) -C $(MANUAL_DIR) breathe_autogen

sphinx-docs-html:
	@$(MAKE) -C $(MANUAL_DIR) sphinx-html

sphinx-docs-latex:
	@$(MAKE) -C $(MANUAL_DIR) sphinx-pdf

manual:
	@$(MAKE) -C $(MANUAL_DIR)

manual-html:
	@$(MAKE) -C $(MANUAL_DIR) html

manual-pdf:
	@$(MAKE) -C $(MANUAL_DIR) pdf

manual-man:
	@$(MAKE) -C $(MANUAL_DIR) man

manual-all:
	@$(MAKE) -C $(MANUAL_DIR) man pdf html
