all: boutcore.so test
PY?=python3
boutcore.pyx: boutcore.pyx.in
setup.py:setup.py.in
resolve_enum.pxd:resolve_enum.pxd.in
helper.h:helper.h.in
helper.cxx:helper.cxx.in
boutcpp.pxd:boutcpp.pxd.in

TOGEN=setup.py boutcore.pyx resolve_enum.pxd helper.cxx helper.h boutcpp.pxd
$(TOGEN): Makefile
	@echo generating $@
	@PATH=../../../bin:$$PATH PY=$(PY) bash $@.in > $@

boutcore.so: boutcore.pyx setup.py helper.cxx helper.h $(TOGEN)
	@$(PY) setup.py build_ext --inplace

test: boutcore.so
	@PYTHONPATH=..:$$PYTHONPATH $(PY) ./runtest > run.log

check: test
.PHONY:

clean: .PHONY
	rm -f boutcore.cpp boutcore.so
	rm -f $(TOGEN)
	rm -rf build

pseudoinstall: check
	cp boutcore*.so ..
