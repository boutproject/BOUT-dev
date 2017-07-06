all: boutcore.so test
ifndef PY
PY=python3
endif
boutcore.pyx: boutcore.pyx.in
setup.py:setup.py.in
resolve_enum.pxd:resolve_enum.pxd.in
helper.h:helper.h.in
helper.cxx:helper.cxx.in

TOGEN=setup.py boutcore.pyx resolve_enum.pxd helper.cxx helper.h
$(TOGEN): Makefile
	@echo generating $@
	@PATH=../../../bin:$(PATH) bash $@.in > $@

boutcore.so: boutcore.pyx setup.py helper.cxx helper.h $(TOGEN)
	@$(PY) setup.py build_ext --inplace

test: boutcore.so
	@$(PY) ./runtest

check: test
.PHONY:

clean: .PHONY
	rm -f boutcore.cpp boutcore.so
	rm -f $(TOGEN)
	rm -rf build

pseudoinstall: check
	cp boutcore*.so ..
