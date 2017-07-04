all: boutcore.so

boutcore.pyx: boutcore.pyx.in resolve_enum.pxd
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
	rm boutcore.cpp
	rm boutcore.so
	rm $(TOGEN)

#FLAGS=-pipe -Wall -DNDEBUG -O2 -g -pipe -Wall -fPIC -fPIC
FLAGS=-O
compile_hand:
	mpic++ -I../../include  -I/usr/include -I/usr/include/hdf -I/usr/include   -g -O2 -O -DCHECK=2 -DSIGHANDLE -DBOUT_VERSION_STRING="\"4.0.0\"" -DBOUT_VERSION_DOUBLE=4.00  -DBACKTRACE -DNCDF -DHDF5 -DLAPACK -DBOUT_HAS_PVODE   -c boutcore.cpp -o build/temp.linux-x86_64-2.7/boutcore.o -I/usr/include/python2.7 
	mpic++ -I../../include  -I/usr/include -I/usr/include/hdf -I/usr/include   -g -O2 -O -DCHECK=2 -DSIGHANDLE -DBOUT_VERSION_STRING="\"4.0.0\"" -DBOUT_VERSION_DOUBLE=4.00  -DBACKTRACE -DNCDF -DHDF5 -DLAPACK -DBOUT_HAS_PVODE   -c helper.cxx -o build/temp.linux-x86_64-2.7/helper.o -I/usr/include/python2.7 
#	mpic++ $(FLAGS) -I../../include -I/usr/include/python2.7 -c boutcore.cpp -o build/temp.linux-x86_64-2.7/boutcore.o
# -std=c++14
#	mpic++ $(FLAGS) -I../../include -I/usr/include/python2.7 -c helper.cxx -o build/temp.linux-x86_64-2.7/helper.o
#-std=c++14
#-pthread -shared -Wl,-z,relro 
	mpic++ -shared  build/temp.linux-x86_64-2.7/boutcore.o build/temp.linux-x86_64-2.7/helper.o -L../../lib -L/usr/lib64 -lm -lbout++ -lfftw3 -lnetcdf_c++ -lnetcdf -lm -ldl -lz -lhdf5_hl -lhdf5 -lhdf5_hl -lblas -llapack -lpvode -lpvpre -lpython2.7 -o /home/dave/soft/BOUT-dev/issue-34/examples/cython/boutcore.so
