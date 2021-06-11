# Makefile for example codes
# Works in the same way as the root Makefile

BOUT_TOP = ../..

BUILD = $(shell PYTHONPATH=../../tools/pylib:$$PYTHONPATH ./test_suite_make --get-list)
BUILD_ALL = $(shell PYTHONPATH=../../tools/pylib:$$PYTHONPATH ./test_suite_make --get-list --all)

CHECKING = $(shell PYTHONPATH=../../tools/pylib:$$PYTHONPATH ./test_suite --get-list)
CHECKING_ALL = $(shell PYTHONPATH=../../tools/pylib:$$PYTHONPATH ./test_suite --get-list --all)

include $(BOUT_TOP)/make.config

CHECKING_ = $(CHECKING:%=%_checking)
CHECKING_ALL_ = $(CHECKING_ALL:%=%_checking)

DIRS_CLEAN = $(BUILD)

.PHONY: $(BUILD_ALL) $(CHECKING_ALL) check

all: $(BUILD)

all-all: $(BUILD_ALL)

buildncheck: all check

$(BUILD):
	@$(MAKE) --no-print-directory -C $@

check: $(CHECKING_)
	@echo $(CHECKING_)

$(CHECKING_ALL_):
	@$(MAKE) --no-print-directory -C $(@:%_checking=%) runtest
