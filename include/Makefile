all: derivs.hxx

derivs.hxx: derivs.hxx.in
	sh $< > $@ || (rm $@ ; exit 1)
