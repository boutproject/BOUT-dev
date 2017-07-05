all: boutcore.so

boutcore.pyx: boutcore.pyx.in
setup.py:setup.py.in
resolve_enum.pxd:resolve_enum.pxd.in
helper.h:helper.h.in
helper.cxx:helper.cxx.in

TOGEN= setup.py boutcore.pyx resolve_enum.pxd helper.cxx helper.h
$(TOGEN): Makefile
	@echo generating $@
	@PATH=../../bin:$(PATH) bash $@.in > $@

boutcore.so: boutcore.pyx setup.py helper.cxx helper.h $(TOGEN)
	@python setup.py build_ext --inplace

test: all
	@python ./runtest

check: test
.PHONY:

clean: .PHONY
	rm boutcore.cpp boutcore.so
	rm $(TOGEN)

