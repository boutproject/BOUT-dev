
BOUT_TOP	= .

DIRS			= src

ifndef TARGET
TARGET=libfast
endif
# Add this to DIRS to have examples compiled
#examples

include make.config

######################################################################
# Tests
######################################################################

check-unit-tests:
	@$(MAKE) --no-print-directory -C tests/unit check

check-mms-tests:
	@cd tests/MMS; ./test_suite

check-integrated-tests:
	@cd tests/integrated; ./test_suite_make
	@cd tests/integrated; PYTHONPATH=${PWD}/../../tools/pylib/:${PYTHONPATH} ./test_suite


check: check-unit-tests check-integrated-tests check-mms-tests
