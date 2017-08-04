
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
	@$(MAKE) --no-print-directory -C tests/MMS check

check-integrated-tests:
	@$(MAKE) --no-print-directory -C tests/integrated check

check: check-unit-tests check-integrated-tests check-mms-tests
