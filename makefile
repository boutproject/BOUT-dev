
BOUT_TOP	= .

DIRS			= src

ifndef TARGET
TARGET=libfast
endif
# Add this to DIRS to have examples compiled
#examples

include make.config

check:
	@$(MAKE) --no-print-directory -C examples check
