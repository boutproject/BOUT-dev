all: field3d.so

field3d.so: field3d.pyx setup.py
	@python setup.py build_ext --inplace

test: all
	@python test.py
