BOUT_TOP = .

include make.config

ifneq ("$(PVODE)foo", "foo")
	# libraries then main source directory
	srcdirs = PVODE/source PVODE/precon
	# manual
endif
srcdirs += SRC


all: PVODE/lib make.config $(srcdirs)

.PHONY:$(srcdirs)
$(srcdirs):
	@$(MAKE) -C $@

config:
	./configure

clean::
	@echo Cleaning...
	@for pp in $(srcdirs); do $(MAKE) -C $$pp clean; done
	@$(RM) include/*
	@$(MAKE) -C examples clean;

PVODE/lib:
	$(MKDIR) $@

.PHONY: examples
examples: all
	@$(MAKE) -C examples;