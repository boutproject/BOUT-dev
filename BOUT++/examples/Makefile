# Makefile for example codes
# Works in the same way as the root Makefile

testdirs = Advect1D RayleighTaylor DriftInstability Uedge_benchmark

.PHONY: all
all: $(testdirs)

.PHONY: $(testdirs)
$(testdirs):
	$(MAKE) -C $@;

.PHONY: clean
clean:
	@for pp in $(testdirs); do $(MAKE) -C $$pp clean; done
