all: boutcore.so
PY?=python3
BOUT_TOP?=../../..
boutcore.pyx: boutcore.pyx.in helper.py resolve_enum_inv.pyx.in
setup.py:setup.py.in
resolve_enum.pxd:resolve_enum.pxd.in
helper.h:helper.h.in
helper.cxx:helper.cxx.in
boutcpp.pxd:boutcpp.pxd.in

TOGEN=setup.py boutcore.pyx resolve_enum.pxd helper.cxx helper.h boutcpp.pxd
$(TOGEN): Makefile $(BOUT_TOP)/make.config $(BOUT_TOP)/bin/bout-config common.sh
	@echo "  Generating $@"
	@PATH=$(BOUT_TOP)/bin:$$PATH PY=$(PY) bash $@.in > $@.tmp \
		|| (fail=$$?; echo "touch $@ to ignore failed generation" ; exit $$fail)
	@mv $@.tmp $@

boutcore.so: boutcore.pyx setup.py helper.cxx helper.h $(TOGEN)
	@$(PY) setup.py build_ext --inplace

.PHONY: clean

clean:
	rm -f boutcore.cpp boutcore.so
	rm -f $(TOGEN)
	rm -rf build

pseudoinstall: boutcore.so
	cp boutcore*.so ..

.PHONY: input
input:
	@PYTHONPATH=..:$$PYTHONPATH $(PY) geninput.py

.PHONY: debug
debug: debug.so
	$(PY) ./debug.py
	! ./debug.py
	! mpirun -np 1 $(PY) ./debug.py
	! mpirun -np 1 ./debug.py

debug.so: debug.pyx debug.pxd debug_setup.py
	$(PY) ./debug_setup.py build_ext --inplace
