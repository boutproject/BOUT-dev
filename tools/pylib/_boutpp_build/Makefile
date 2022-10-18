all: boutpp.so
PY?=python3
BOUT_TOP?=../../..

boutpp.pyx: boutpp.pyx.jinja helper.py resolve_enum_inv.pyx.jinja scan_enums.py
resolve_enum.pxd:resolve_enum.pxd.jinja scan_enums.py
helper.h:helper.h.jinja
helper.cxx:helper.cxx.jinja
boutcpp.pxd:boutcpp.pxd.jinja

TOGEN=helper.cxx helper.h boutcpp.pxd boutpp.pyx resolve_enum.pxd
$(TOGEN): Makefile $(BOUT_TOP)/make.config $(BOUT_TOP)/bin/bout-config ../boutconfig/__init__.py

$(TOGEN): generate.py
	@echo "  Generating $@"
	@PYTHONPATH=..:$$PYTHONPATH $(PY) generate.py $@.jinja $@

# Setting LDSHARED sets the linker. However, setting it might drop
# flags - thus it is not set by default. If mpic++ -std=c++14 is
# different to the python compiler, LDSHARED might need to be set. In
# that case it might be required to check the linking commands, and
# add the apropriate flags.
# LDSHARED=$$($(BOUT_TOP)/bin/bout-config --ld )
boutpp.so: $(TOGEN) $(TOGEN2)
	CC=$$($(BOUT_TOP)/bin/bout-config --cc ) \
	 CXX=$$($(BOUT_TOP)/bin/bout-config --cxx) \
	 PYTHONPATH=..:$$PYTHONPATH $(PY) setup.py build_ext --inplace

.PHONY: clean

clean:
	rm -f boutpp.cpp boutpp.so
	rm -f $(TOGEN) $(TOGEN2)
	rm -rf build

pseudoinstall: boutpp.so
	cp boutpp*.so ..

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
