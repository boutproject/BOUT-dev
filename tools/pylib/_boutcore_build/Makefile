all: boutcore.so test
PY?=python3
boutcore.pyx: boutcore.pyx.in helper.py resolve_enum_inv.pyx.in
setup.py:setup.py.in
resolve_enum.pxd:resolve_enum.pxd.in
helper.h:helper.h.in
helper.cxx:helper.cxx.in
boutcpp.pxd:boutcpp.pxd.in

TOGEN=setup.py boutcore.pyx resolve_enum.pxd helper.cxx helper.h boutcpp.pxd
$(TOGEN): Makefile
	@echo "  Generating $@"
	@PATH=../../../bin:$$PATH PY=$(PY) bash $@.in > $@.tmp \
		|| (fail=$?; echo "touch $@ to ignore failed generation" ; exit $fail)
	@mv $@.tmp $@

boutcore.so: boutcore.pyx setup.py helper.cxx helper.h $(TOGEN)
	@$(PY) setup.py build_ext --inplace

test: boutcore.so
	@PYTHONPATH=..:$$PYTHONPATH $(PY) ./runtest 2>&1 > run.log || (ex=$$? ; cat run.log ; exit $$ex)
	@PYTHONPATH=..:$$PYTHONPATH $(PY) mini.py 2>&1 > run.log || (ex=$$? ; cat run.log ; exit $$ex)

check: test
.PHONY: clean

clean:
	rm -f boutcore.cpp boutcore.so
	rm -f $(TOGEN)
	rm -rf build

pseudoinstall: check
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

debug.so: debug.pyx debug.pxd setup_debug.py
	$(PY) ./setup_debug.py build_ext --inplace
