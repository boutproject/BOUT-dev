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
	@$(MAKE) --no-print-directory -C $@

config:
	./configure

clean::
	@for pp in $(srcdirs); do echo Cleaning $$pp; $(MAKE) --no-print-directory -C $$pp clean; done
	@echo include cleaned
	@$(RM) include/*
	@echo lib cleaned
	@$(RM) lib/*
	@echo Cleaning examples
	@$(MAKE) --no-print-directory -C examples clean;

PVODE/lib:
	$(MKDIR) $@

.PHONY: examples
examples: all
	@$(MAKE) -C examples;
