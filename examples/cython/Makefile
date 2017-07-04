all: field3d.so

field3d.so: field3d.pyx setup.py helper.cxx helper.h
	@python setup.py build_ext --inplace

test: all
	@python test.py

check: test
.PHONY:

clean: .PHONY
	rm field3d.cpp
	rm field3d.so

#FLAGS=-pipe -Wall -DNDEBUG -O2 -g -pipe -Wall -fPIC -fPIC
FLAGS=-O
compile_hand:
	mpic++ -I../../include  -I/usr/include -I/usr/include/hdf -I/usr/include   -g -O2 -O -DCHECK=2 -DSIGHANDLE -DBOUT_VERSION_STRING="\"4.0.0\"" -DBOUT_VERSION_DOUBLE=4.00  -DBACKTRACE -DNCDF -DHDF5 -DLAPACK -DBOUT_HAS_PVODE   -c field3d.cpp -o build/temp.linux-x86_64-2.7/field3d.o -I/usr/include/python2.7 
	mpic++ -I../../include  -I/usr/include -I/usr/include/hdf -I/usr/include   -g -O2 -O -DCHECK=2 -DSIGHANDLE -DBOUT_VERSION_STRING="\"4.0.0\"" -DBOUT_VERSION_DOUBLE=4.00  -DBACKTRACE -DNCDF -DHDF5 -DLAPACK -DBOUT_HAS_PVODE   -c helper.cxx -o build/temp.linux-x86_64-2.7/helper.o -I/usr/include/python2.7 
#	mpic++ $(FLAGS) -I../../include -I/usr/include/python2.7 -c field3d.cpp -o build/temp.linux-x86_64-2.7/field3d.o
# -std=c++14
#	mpic++ $(FLAGS) -I../../include -I/usr/include/python2.7 -c helper.cxx -o build/temp.linux-x86_64-2.7/helper.o
#-std=c++14
#-pthread -shared -Wl,-z,relro 
	mpic++ -shared  build/temp.linux-x86_64-2.7/field3d.o build/temp.linux-x86_64-2.7/helper.o -L../../lib -L/usr/lib64 -lm -lbout++ -lfftw3 -lnetcdf_c++ -lnetcdf -lm -ldl -lz -lhdf5_hl -lhdf5 -lhdf5_hl -lblas -llapack -lpvode -lpvpre -lpython2.7 -o /home/dave/soft/BOUT-dev/issue-34/examples/cython/field3d.so
