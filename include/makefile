all: derivs.hxx


derivs.hxx: derivs.hxx.in.py derivs.hxx.in.function derivs.hxx.in.description
# run first dependency and write output to target
# on failure delete output and report error
	./$< > $@ || (rm $@ ; exit 1)
