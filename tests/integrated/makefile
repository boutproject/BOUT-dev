# Makefile for example codes
# Works in the same way as the root Makefile

BOUT_TOP = ../..

BUILD = $(shell PYTHONPATH=../../tools/pylib:$$PYTHONPATH ./test_suite_make --get-list)
BUILDALL = $(shell PYTHONPATH=../../tools/pylib:$$PYTHONPATH ./test_suite_make --get-list --all)

CHECKING = $(shell PYTHONPATH=../../tools/pylib:$$PYTHONPATH ./test_suite --get-list)

include $(BOUT_TOP)/make.config

CHECKING_ = $(CHECKING:%=%_checking)

DIRS_CLEAN = $(BUILD)

.PHONY: $(BUILDALL) $(CHECKING_) check

all: $(BUILD)

allall: $(BUILDALL)

buildncheck: all check

$(BUILDALL):
	@$(MAKE) --no-print-directory -C $@

check: $(CHECKING_)

$(CHECKING_):
	@$(MAKE) --no-print-directory -C $(@:%_checking=%) runtest
